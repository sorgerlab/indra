"""Read a list of files located in your local directory."""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import logging
import pickle
import random

logger = logging.getLogger('indra.tools.reading.read_files')

from indra.tools.reading.util.script_tools import get_parser
from indra.tools.reading.readers import get_dir, get_reader_classes, Content


def make_parser():
    """Create the argument parser, derived from the general scripts parser."""
    parser = get_parser(
        __doc__,
        ('A file containing a list of files/file paths to be read. These '
         'should be nxml or txt files.')
    )
    parser.add_argument(
        dest='output_name',
        help=('Results will be pickled in files '
              '<output_name>_stmts.pkl and <output_name>_readings.pkl.')
    )
    return parser


def read_files(files, readers, **kwargs):
    """Read the files in `files` with the reader objects in `readers`.

    Parameters
    ----------
    files : list [str]
        A list of file paths to be read by the readers. Supported files are
        limited to text and nxml files.
    readers : list [Reader instances]
        A list of Reader objects to be used reading the files.
    **kwargs :
        Other keyword arguments are passed to the `read` method of the readers.

    Returns
    -------
    output_list : list [ReadingData]
        A list of ReadingData objects with the contents of the readings.
    """
    reading_content = [Content.from_file(filepath) for filepath in files]
    output_list = []
    for reader in readers:
        res_list = reader.read(reading_content, **kwargs)
        if res_list is None:
            logger.info("Nothing read by %s." % reader.name)
        else:
            logger.info("Successfully read %d content entries with %s."
                        % (len(res_list), reader.name))
            output_list += res_list
    logger.info("Read %s text content entries in all." % len(output_list))
    return output_list


def main():
    # Load arguments.
    parser = make_parser()
    args = parser.parse_args()
    if args.debug and not args.quiet:
        logger.setLevel(logging.DEBUG)

    # Load the input file.
    with open(args.input_file, 'r') as f:
        input_lines = f.readlines()
    logger.info("Found %d files." % len(input_lines))
    for ftype in ['nxml', 'txt']:
        logger.debug('%d are %s' % (
            len([f for f in input_lines if f.endswith(ftype)]), ftype
        ))

    # Select only a sample of the lines, if sample is chosen.
    if args.n_samp is not None:
        input_lines = random.sample(input_lines, args.n_samp)

    # If a range is specified, only use that range.
    if args.range_str is not None:
        start_idx, end_idx = [int(n) for n in args.range_str.split(':')]
        input_lines = input_lines[start_idx:end_idx]

    # Create a single base directory
    base_dir = get_dir('run_%s' % ('_and_'.join(args.readers)))

    # Set the verbosity. The quiet argument overrides the verbose argument.
    verbose = args.verbose and not args.quiet

    # Get the readers objects.
    readers = [reader_class(base_dir=base_dir, n_proc=args.n_proc)
               for reader_class in get_reader_classes()
               if reader_class.name.lower() in args.readers]

    # Read the files.
    outputs = read_files(input_lines, readers, verboes=verbose)
    reading_out_path = args.name + '_readings.pkl'
    with open(reading_out_path, 'wb') as f:
        pickle.dump([output.make_tuple(None) for output in outputs], f)
    print("Reading outputs stored in %s." % reading_out_path)

    stmts = [s for rd in outputs for s in rd.get_statements()]
    stmts_pkl_path = args.name + '_stmts.pkl'
    with open(stmts_pkl_path, 'wb') as f:
        pickle.dump(stmts, f)
        print("Statements pickled in %s." % stmts_pkl_path)


if __name__ == '__main__':
    main()
