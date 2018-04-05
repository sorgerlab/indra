from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, int

import logging
from os import path
from functools import wraps
from datetime import datetime

logger = logging.getLogger('reading_manager')
THIS_DIR = path.dirname(path.abspath(__file__))


class ReadingManager(object):
    reader = NotImplemented

    def __init__(self, buffer_days=1):
        self.buffer_days = buffer_days
        self.reader_version = None
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
                is_init_upload = (func.__name__ == 'populate')
                db.insert('reading_updates', init_upload=is_init_upload,
                          source=self.my_source,
                          run_datetime=self.run_datetime,
                          earliest_datetime=self.begin_datetime,
                          latest_datetime=self.end_datetime)
            return completed
        return run_and_record_update

    def _get_latest_update(self, db):
        """Get the date of the latest update."""
        update_list = db.select_all(db.ReadingUpdates,
                                    db.ReadingUpdates.reader == self.reader)
        if not len(update_list):
            logger.error("The database has not had an initial upload, or else "
                         "the updates table has not been populated.")
            return False

        return max([u.datetime for u in update_list])

    @ReadingManager._handle_update_table
    def read_all(self, db, n_proc):
        """Read everything available on the database."""

    @ReadingManager._handle_update_table
    def read_new(self, db, n_proc):
        """Update the readings and raw statements in the database."""


class ReachManager(ReadingManager):
    reader = 'reach'


class SparserManager(ReadingManager):
    reader = 'sparser'
