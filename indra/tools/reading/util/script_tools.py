"""Useful tools for reading scripts."""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import logging
import json
from argparse import ArgumentParser
from indra.util.get_version import get_version as get_indra_version
from multiprocessing import Pool

logger = logging.getLogger('script_tools')


def get_parser(description, input_desc):
    """Get a parser that is generic to reading scripts.

    Parameters
    ----------
    description : str
        A description of the tool, usually about one line long.
    input_desc: str
        A string describing the nature of the input file used by the reading
        tool.

    Returns
    -------
    parser : argparse.ArgumentParser instance
        An argument parser object, to which further arguments can be added.
    """
    parser = ArgumentParser(description=description)
    parser.add_argument(
        dest='input_file',
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
        dest='n_proc',
        help='Select the number of processes to use.',
        type=int,
        default=1
        )
    parser.add_argument(
        '-s', '--sample',
        dest='n_samp',
        help='Read a random sample of size N_SAMP of the inputs.',
        type=int
        )
    parser.add_argument(
        '-I', '--in_range',
        dest='range_str',
        help='Only read input lines in the range given as <start>:<end>.'
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
    # parser.add_argument(
    #     '-m', '--messy',
    #     help='Do not clean up directories created while reading.',
    #     action='store_true'
    #     )
    return parser


class StatementData(object):
    """Contains metadata for statements, as well as the statement itself.

    This, like ReadingData, is primarily designed for use with the database,
    carrying valuable information and methods for such.

    Parameters
    ----------
    statement : an indra Statement instance
        The statement whose extra meta data this object encapsulates.
    reading_id : int or None
        The id number of the entry in the `readings` table of the database.
        None if no such id is available.
    """
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


def get_stmts_safely(reading_data):
    stmt_data_list = []
    try:
        stmts = reading_data.get_statements()
    except Exception as e:
        logger.error("Got exception creating statements for %d."
                     % reading_data.reading_id)
        logger.exception(e)
        return
    if stmts is not None:
        for stmt in stmts:
            stmt_data_list.append(StatementData(stmt, reading_data.reading_id))
    return stmt_data_list


def make_statements(reading_data_list, num_proc=1):
    """Convert a list of ReadingData instances into StatementData instances."""
    stmt_data_list = []

    if num_proc is 1:  # Don't use pool if not needed.
        for reading_data in reading_data_list:
            stmt_data_list += get_stmts_safely(reading_data)
    else:
        try:
            pool = Pool(num_proc)
            stmt_data_list_list = pool.map(get_stmts_safely, reading_data_list)
            for stmt_data_sublist in stmt_data_list_list:
                if stmt_data_sublist is not None:
                    stmt_data_list += stmt_data_sublist
        finally:
            pool.close()
            pool.join()

    logger.info("Found %d statements from %d readings." %
                (len(stmt_data_list), len(reading_data_list)))
    return stmt_data_list
