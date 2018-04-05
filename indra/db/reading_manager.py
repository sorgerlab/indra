from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, int

import logging
from os import path
from functools import wraps
from datetime import datetime

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
        self.buffer_days = buffer_days
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

    def _get_latest_update(self, db):
        """Get the date of the latest update."""
        update_list = db.select_all(
            db.ReadingUpdates,
            db.ReadingUpdates.reader == self.reader.name
            )
        if not len(update_list):
            logger.error("The database has not had an initial upload, or else "
                         "the updates table has not been populated.")
            return False

        return max([u.datetime for u in update_list])


class BulkReadingManager(ReadingManager):
    @ReadingManager._handle_update_table
    def read_all(self, db, n_proc):
        """Read everything available on the database."""
        self.end_datetime = self.run_datetime
        return True

    @ReadingManager._handle_update_table
    def read_new(self, db, n_proc):
        """Update the readings and raw statements in the database."""
        return True
