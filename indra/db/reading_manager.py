from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, int

import logging
from os import path
from functools import wraps
from datetime import datetime, timedelta
from indra.db.util import get_primary_db, get_test_db

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
    def run_reading(self, db, id_dict, n_proc, verbose):
        base_dir = path.join(THIS_DIR, 'read_all_%s' % self.reader.name)
        reader_inst = self.reader(base_dir=base_dir, n_proc=n_proc)

        logger.info("Making readings...")
        outputs = rdb.produce_readings(id_dict, [reader_inst], verbose=verbose,
                                       read_mode='unread_unread', db=db,
                                       prioritize=True)
        logger.info("Made %d readings." % len(outputs))
        logger.info("Making statements...")
        rdb.produce_statements(outputs, n_proc=n_proc, db=db)
        return

    @ReadingManager._handle_update_table
    def read_all(self, db, n_proc=1, verbose=True):
        """Read everything available on the database."""
        self.end_datetime = self.run_datetime
        trids = {trid for trid, in db.select_all(db.TextContent.text_ref_id)}
        logger.info("Producing readings for %d text refs." % len(trids))
        id_dict = {'trid': trids}
        self.run_reading(db, id_dict, n_proc, verbose)
        return True

    @ReadingManager._handle_update_table
    def read_new(self, db, n_proc=1, verbose=True):
        """Update the readings and raw statements in the database."""
        self.end_datetime = self.run_datetime
        self.begin_datetime = self._get_latest_updatetime(db) - self.buffer
        trid_q = db.filter_query(
            db.TextContent.text_ref_id,
            db.TextContent.insert_date > self.begin_datetime
            )
        id_dict = {'trid': {trid for trid, in trid_q.all()}}
        logger.info("Producing readings for %d new text refs."
                    % len(id_dict['trid']))
        self.run_reading(db, id_dict, n_proc, verbose)
        return True


if __name__ == '__main__':
    if args.test:
        db = get_test_db()
    else:
        db = get_primary_db()

    bulk_managers = [BulkReadingManager(reader_name, buffer_days=args.buffer)
                     for reader_name in ['REACH', 'SPARSER']]

    if args.task == 'read_all':
        for bulk_manager in bulk_managers:
            bulk_manager.read_all(db, n_proc=args.num_procs)
    elif args.task == 'read_new':
        for bulk_manager in bulk_managers:
            bulk_manager.read_new(db, n_proc=args.num_procs)
