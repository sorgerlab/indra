from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, int

import logging
from os import path
from functools import wraps
from datetime import datetime, timedelta
from indra.db.util import get_primary_db, get_test_db, get_db
from indra.tools.reading.submit_reading_pipeline import submit_db_reading,\
    wait_for_complete

if __name__ == '__main__':
    # NOTE: PEP8 will complain about this, however having the args parsed up
    # here prevents a long wait just to fined out you entered a command wrong.
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(
        description='Manage content on INDRA\'s database.'
        )
    parent_read_parser = ArgumentParser(add_help=False)
    parent_read_parser.add_argument(
        choices=['read_all', 'read_new'],
        dest='task',
        help=('Choose whether you want to read/reread everything, or only '
              'read the content added since the last update. Note that '
              'content from one day before the latests update will also be '
              'checked, to avoid content update overlap errors.')
        )
    parent_read_parser.add_argument(
        '-t', '--test',
        action='store_true',
        help='Run tests using one of the designated test databases.'
        )
    parent_read_parser.add_argument(
        '-b', '--buffer',
        type=int,
        default=1,
        help=('Set the number number of buffer days read prior to the most '
              'recent update. The default is 1 day.')
        )
    local_read_parser = ArgumentParser(add_help=False)
    local_read_parser.add_argument(
        '-n', '--num_procs',
        dest='num_procs',
        type=int,
        default=1,
        help=('Select the number of processors to use during this operation. '
              'Default is 1.')
        )
    parser.add_argument(
        '--use_batch',
        action='store_true',
        help=('Choose to run the update on amazon batch. Note that this '
              'option will cause the -n/--num_procs option to be ignored.')
        )
    parser.add_argument(
        '--database',
        default='primary',
        help=('Select a database from the names given in the config or '
              'environment, for example primary is INDRA_DB_PRIMAY in the '
              'config file and INDRADBPRIMARY in the environment. The default '
              'is \'primary\'. Note that this is overwridden by use of the '
              '--test flag if \'test\' is not a part of the name given.')
        )
    aws_read_parser = ArgumentParser(add_help=False)
    aws_read_parser.add_argument(
        '--project_name',
        help=('For use with --use_batch. Set the name of the project for '
              'which this reading is being done. This is used to label jobs '
              'on aws batch for monitoring and accoutning purposes.')
        )
    subparsers = parser.add_subparsers(title='Method')
    subparsers.required = True
    subparsers.dest = 'method'

    local_desc = 'Run the reading for the update locally.'
    subparsers.add_parser(
        'local',
        parents=[parent_read_parser, local_read_parser],
        help=local_desc,
        description=local_desc,
        formatter_class=ArgumentDefaultsHelpFormatter
        )
    aws_desc = 'Run the reading for the update on amazon batch.'
    subparsers.add_parser(
        'aws',
        parents=[parent_read_parser, aws_read_parser],
        help=local_desc,
        description=local_desc,
        formatter_class=ArgumentDefaultsHelpFormatter
        )
    args = parser.parse_args()

from indra.tools.reading.db_reading import read_db as rdb
from indra.tools.reading.readers import get_reader_class

logger = logging.getLogger('reading_manager')
THIS_DIR = path.dirname(path.abspath(__file__))


class ReadingUpdateError(Exception):
    pass


class ReadingManager(object):
    """Abstract class for managing the readings of the database.

    Parameters
    ----------
    reader_name : str
        The name of the reader to be used in a given run of reading.
    buffer_days : int
        The number of days before the previous update/initial upload to look for
        "new" content to be read. This prevents any issues with overlaps between
        the content upload pipeline and the reading pipeline.
    """
    def __init__(self, reader_name, buffer_days=1):
        self.reader = get_reader_class(reader_name)
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
            logger.warning("The database has not had an initial upload, or "
                           "else the updates table has not been populated.")
            return None

        return max([u.latest_datetime for u in update_list])

    def read_all(self, db):
        """Perform an initial reading all content in the database (populate).

        This must be defined in a child class.
        """
        raise NotImplementedError

    def read_new(self, db):
        """Read only new content (update).

        This must be defined in a child class.
        """
        raise NotImplementedError


class BulkReadingManager(ReadingManager):
    """An abstract class which defines methods required for reading in bulk.

    This takes exactly the parameters used by :py:class:`ReadingManager`.
    """
    def _run_reading(self, db, trids, max_refs=5000):
        raise NotImplementedError("_run_reading must be defined in child.")

    @ReadingManager._handle_update_table
    def read_all(self, db):
        """Read everything available on the database."""
        self.end_datetime = self.run_datetime
        trids = {trid for trid, in db.select_all(db.TextContent.text_ref_id)}
        self._run_reading(db, trids)
        return True

    @ReadingManager._handle_update_table
    def read_new(self, db):
        """Update the readings and raw statements in the database."""
        self.end_datetime = self.run_datetime
        latest_updatetime = self._get_latest_updatetime(db)
        if latest_updatetime is not None:
            self.begin_datetime = latest_updatetime - self.buffer
        else:
            raise ReadingUpdateError("There are no previous updates. "
                                     "Please run_all.")
        trid_q = db.filter_query(
            db.TextContent.text_ref_id,
            db.TextContent.insert_date > self.begin_datetime
            )
        trids = {trid for trid, in trid_q.all()}
        self._run_reading(db, trids)
        return True


class BulkAwsReadingManager(BulkReadingManager):
    """This is the reading manager when updating using AWS Batch.

    This takes all the parameters used by :py:class:`BulkReadingManager`, and
    in addition:

    Parameters
    ----------
    project_name : str
        You can select a name for the project for which this reading is being
        run. This name has a default value set in your config file. The batch
        jobs used in reading will be tagged with this project name, for
        accounting purposes.
    """
    def __init__(self, *args, **kwargs):
        self.project_name = kwargs.pop('project_name', None)
        super(BulkAwsReadingManager, self).__init__(*args, **kwargs)
        return

    def _run_reading(self, db, trids, max_refs=5000):
        if len(trids)/max_refs >= 1000:
            raise ReadingUpdateError("Too many id's for one submission. "
                                     "Break it up and do it manually.")

        logger.info("Producing readings on aws for %d text refs with new "
                    "content not read by %s." % (len(trids), self.reader.name))
        job_prefix = ('%s_reading_%s'
                      % (self.reader.name.lower(),
                         self.run_datetime.strftime('%Y%m%d_%H%M%S')))
        with open(job_prefix + '.txt', 'w') as f:
            f.write('\n'.join(['trid:%s' % trid for trid in trids]))
        logger.info("Submitting jobs...")
        job_ids = submit_db_reading(job_prefix, job_prefix + '.txt',
                                    readers=[self.reader.name.lower()],
                                    start_ix=0, end_ix=None,
                                    pmids_per_job=max_refs, num_tries=2,
                                    force_read=False, force_fulltext=False,
                                    read_all_fulltext=False,
                                    project_name=self.project_name)
        logger.info("Waiting for complete...")
        wait_for_complete('run_db_reading_queue', job_list=job_ids,
                          job_name_prefix=job_prefix,
                          idle_log_timeout=1200,
                          kill_on_log_timeout=True,
                          stash_log_method='s3')
        return


class BulkLocalReadingManager(BulkReadingManager):
    """This is the reading manager to be used when running reading locally.

    This takes all the parameters used by :py:class:`BulkReadingManager`, and
    in addition:

    Parameters
    ----------
    n_proc : int
        The number of processed to dedicate to reading. Note the some of the
        readers (e.g. REACH) do not always obey these restrictions.
    verbose : bool
        If True, more detailed logs will be printed. Default is False.
    """
    def __init__(self, *args, **kwargs):
        self.n_proc = kwargs.pop('n_proc', 1)
        self.verbose = kwargs.pop('verbose', False)
        super(BulkLocalReadingManager, self).__init__(*args, **kwargs)
        return

    def _run_reading(self, db, trids, max_refs=5000):
        if len(trids) > max_refs:
            raise ReadingUpdateError("Too many id's to run locally. Try "
                                     "running on batch (use_batch).")
        logger.info("Producing readings locally for %d new text refs."
                    % len(trids))
        base_dir = path.join(THIS_DIR, 'read_all_%s' % self.reader.name)
        reader_inst = self.reader(base_dir=base_dir, n_proc=self.n_proc)

        logger.info("Making readings...")
        outputs = rdb.produce_readings({'trid': trids}, [reader_inst],
                                       read_mode='unread_unread', db=db,
                                       prioritize=True, verbose=self.verbose)
        logger.info("Made %d readings." % len(outputs))
        logger.info("Making statements...")
        rdb.produce_statements(outputs, n_proc=self.n_proc, db=db)
        return


if __name__ == '__main__':
    if args.test:
        if 'test' not in args.database:
            db = get_test_db()
        else:
            db = get_db(args.database)
    elif args.database == 'primary':
        db = get_primary_db()
    else:
        db = get_db(args.database)

    if args.method == 'local':
        bulk_managers = [BulkLocalReadingManager(reader_name,
                                                 buffer_days=args.buffer,
                                                 n_proc=args.num_procs)
                         for reader_name in ['SPARSER', 'REACH']]
    elif args.method == 'aws':
        bulk_managers = [BulkAwsReadingManager(reader_name,
                                               buffer_days=args.buffer,
                                               project_name=args.project_name)
                         for reader_name in ['SPARSER', 'REACH']]

    for bulk_manager in bulk_managers:
        if args.task == 'read_all':
            bulk_manager.read_all(db)
        elif args.task == 'read_new':
            bulk_manager.read_new(db)
