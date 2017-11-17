"""Useful tools for reading scripts."""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import logging
import json
from argparse import ArgumentParser
from indra.util.get_version import get_version as get_indra_version

logger = logging.getLogger('script_tools')


def get_parser(description, input_desc):
    """Get a parser that is generic to reading scripts."""
    parser = ArgumentParser(description=description)
    parser.add_argument(
        '-i', '--input_file',
        help=input_desc
        )
    parser.add_argument(
        '-r', '--readers',
        choices=['reach', 'sparser'],
        help='List of readers to be used.',
        nargs='+'
        )
    parser.add_argument(
        '-n', '--num_procs',
        help='Select the number of processes to use.',
        type=int,
        default=1
        )
    parser.add_argument(
        '-s', '--sample',
        help='To read just a sample of the pmids in the file.',
        type=int
        )
    parser.add_argument(
        '-I', '--in_range',
        help='A range of pmids to be read in the form <start>:<end>.'
        )
    parser.add_argument(
        '-v', '--verbose',
        help='Include output from the readers.',
        action='store_true'
        )
    parser.add_argument(
        '-q', '--quiet',
        help='Suppress most output. Overrides -v and -d options.',
        action='store_true'
        )
    parser.add_argument(
        '-d', '--debug',
        help='Set the logging to debug level.',
        action='store_true'
        )
    parser.add_argument(
        '-m', '--messy',
        help='Do not clean up directories created while reading.',
        action='store_true'
        )
    return parser


class StatementData(object):
    """Contains metadata for statements, as well as the statement itself."""

    def __init__(self, statement, reading_id):
        self.reading_id = reading_id
        self.statement = statement
        self.indra_version = get_indra_version()
        return

    @classmethod
    def get_cols(self):
        """Get the columns for the tuple returned by `make_tuple`."""
        return ('reader_ref', 'uuid', 'type', 'json', 'indra_version')

    def make_tuple(self):
        """Make a tuple for copying into the database."""
        return (self.reading_id, self.statement.uuid,
                self.statement.__class__.__name__,
                json.dumps(self.statement.to_json()),
                self.indra_version)


def make_statements(reading_data_list):
    stmt_data_list = [StatementData(stmt, output.reading_id)
                      for output in reading_data_list
                      for stmt in output.get_statements()]
    logger.info("Found %d statements from %d readings." %
                (len(stmt_data_list), len(reading_data_list)))
    return stmt_data_list
