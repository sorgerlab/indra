"""Useful tools for reading scripts."""
import logging
from argparse import ArgumentParser

from indra.tools.reading.readers import get_reader_classes

logger = logging.getLogger(__name__)


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
        choices=[rc.name.lower() for rc in get_reader_classes()],
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


