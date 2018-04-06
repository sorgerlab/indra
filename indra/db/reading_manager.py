from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, int

import logging
from os import path
from functools import wraps
from datetime import datetime, timedelta
from indra.db.util import get_primary_db, get_test_db
from indra.tools.reading.submit_reading_pipeline import submit_db_reading,\
    wait_for_complete

if __name__ == '__main__':
    # NOTE: PEP8 will complain about this, however having the args parsed up
    # here prevents a long wait just to fined out you entered a command wrong.
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description='Manage content on INDRA\'s database.'
        )
    parser.add_argument(
        choices=['read_all', 'read_new'],
        dest='task',
        help=('Choose whether you want to read/reread everything, or only '
              'read the content added since the last update. Note that content '
              'from one day before the latests update will also be checked, to '
              'avoid content update overlap errors.')
        )
    parser.add_argument(
        '-n', '--num_procs',
        dest='num_procs',
        type=int,
        default=1,
        help=('Select the number of processors to use during this operation. '
              'Default is 1.')
        )
    parser.add_argument(
        '-t', '--test',
        action='store_true',
        help='Run tests using one of the designated test databases.'
        )
    parser.add_argument(
        '-b', '--buffer',
        type=int,
        default=1,
        help=('Set the number number of buffer days read prior to the most '
              'recent update. The default is 1 day.')
        )
    parser.add_argument(
        '--use_batch',
        action='store_true',
        help=('Choose to run the update on amazon batch. Note that this '
              'option will cause the -n/--num_procs option to be ignored.')
        )
    args = parser.parse_args()

from indra.tools.reading.db_reading import read_db as rdb
from indra.tools.reading.readers import get_reader

logger = logging.getLogger('reading_manager')
THIS_DIR = path.dirname(path.abspath(__file__))


class ReadingUpdateError(Exception):
    pass


class ReadingManager(object):
    def __init__(self, reader_name, buffer_days=1):
        self.reader = get_reader(reader_name)
        if self.reader is None:
            raise ReadingUpdateError('Name of reader was not matched to an '
                                     'available reader.')
        self.buffer = timedelta(days=buffer_days)
        self.reader_version = self.reader.get_version()
        self.run_datetime = None
        self.begin_datetime = None
        self.end_datetime = None
        return

    @classmethod
    def _handle_update_table(cls, func):
        @wraps(func)
        def run_and_record_update(self, db, *args, **kwargs):
            self.run_datetime = datetime.utcnow()
            completed = func(self, db, *args, **kwargs)
            if completed:
                is_complete_read = (func.__name__ == 'read_all')
                db.insert('reading_updates', complete_read=is_complete_read,
                          reader=self.reader.name,
                          reader_version=self.reader_version,
                          run_datetime=self.run_datetime,
                          earliest_datetime=self.begin_datetime,
                          latest_datetime=self.end_datetime)
            return completed
        return run_and_record_update

    def _get_latest_updatetime(self, db):
        """Get the date of the latest update."""
        update_list = db.select_all(
            db.ReadingUpdates,
            db.ReadingUpdates.reader == self.reader.name
            )
        if not len(update_list):
            logger.error("The database has not had an initial upload, or else "
                         "the updates table has not been populated.")
            return False

        return max([u.latest_datetime for u in update_list])


class BulkReadingManager(ReadingManager):
    def run_reading(self, db, trids, n_proc, verbose, max_refs=5000,
                    use_batch=False):
        if use_batch:
            if len(trids)/max_refs >= 1000:
                raise ReadingUpdateError("Too many id's for one submission. "
                                         "Break it up and do it manually.")

            logger.info("Producing readings on aws for %d new text refs."
                        % len(trids))
            job_prefix = ('%s_reading_%s'
                          % (self.reader.name.lower(),
                             self.run_datetime.strftime('%Y%m%d_%H%M%S')))
            with open(job_prefix + '.txt', 'w') as f:
                f.write('\n'.join(['trid:%s' % trid for trid in trids]))
            logger.info("Submitting jobs...")
            job_ids = submit_db_reading(job_prefix, job_prefix + '.txt',
                                        [self.reader.name.lower()], 0, None,
                                        max_refs, 2, False, False, False)
            logger.info("Waiting for complete...")
            wait_for_complete('run_db_reading_queue', job_list=job_ids,
                              job_name_prefix=job_prefix,
                              idle_log_timeout=1200,
                              kill_on_log_timeout=True,
                              stash_log_method='s3')
        else:
            if len(trids) > max_refs:
                raise ReadingUpdateError("Too many id's to run locally. Try "
                                         "running on batch (use_batch).")
            logger.info("Producing readings locally for %d new text refs."
                        % len(trids))
            base_dir = path.join(THIS_DIR, 'read_all_%s' % self.reader.name)
            reader_inst = self.reader(base_dir=base_dir, n_proc=n_proc)

            logger.info("Making readings...")
            outputs = rdb.produce_readings({'trid': trids}, [reader_inst],
                                           read_mode='unread_unread', db=db,
                                           prioritize=True, verbose=verbose)
            logger.info("Made %d readings." % len(outputs))
            logger.info("Making statements...")
            rdb.produce_statements(outputs, n_proc=n_proc, db=db)
        return

    @ReadingManager._handle_update_table
    def read_all(self, db, n_proc=1, verbose=True, use_batch=False):
        """Read everything available on the database."""
        self.end_datetime = self.run_datetime
        trids = {trid for trid, in db.select_all(db.TextContent.text_ref_id)}
        self.run_reading(db, trids, n_proc, verbose, use_batch=use_batch)
        return True

    @ReadingManager._handle_update_table
    def read_new(self, db, n_proc=1, verbose=True, use_batch=False):
        """Update the readings and raw statements in the database."""
        self.end_datetime = self.run_datetime
        self.begin_datetime = self._get_latest_updatetime(db) - self.buffer
        trid_q = db.filter_query(
            db.TextContent.text_ref_id,
            db.TextContent.insert_date > self.begin_datetime
            )
        trids = {trid for trid, in trid_q.all()}
        self.run_reading(db, trids, n_proc, verbose, use_batch=use_batch)
        return True


if __name__ == '__main__':
    if args.test:
        db = get_test_db()
    else:
        db = get_primary_db()

    bulk_managers = [BulkReadingManager(reader_name, buffer_days=args.buffer)
                     for reader_name in ['SPARSER', 'REACH']]

    if args.task == 'read_all':
        for bulk_manager in bulk_managers:
            bulk_manager.read_all(db, n_proc=args.num_procs,
                                  use_batch=args.use_batch)
    elif args.task == 'read_new':
        for bulk_manager in bulk_managers:
            bulk_manager.read_new(db, n_proc=args.num_procs,
                                  use_batch=args.use_batch)
