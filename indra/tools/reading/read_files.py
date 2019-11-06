"""Read a list of files located in your local directory."""
import json
import pickle
import random
import logging
from os import path, listdir

from indra.tools.reading.readers.core import dump_readings
from indra.tools.reading.util.script_tools import get_parser
from indra.tools.reading.readers import get_dir, get_reader_classes, Content

logger = logging.getLogger(__name__)


def make_parser():
    """Create the argument parser, derived from the general scripts parser."""
    parser = get_parser(
        __doc__,
        ('Either a file containing a list of files/file paths to be read or '
         'else a directory containing the files to be read. These '
         'should be nxml or txt files. The basenames of the files will be '
         'used as the IDs for the content.')
    )
    parser.add_argument(
        dest='output_path',
        help=('The location of the results. Results will be stored in files '
              '<output_path>/stmts.json and <output_path>/readings.json, '
              'or likewise with .pkl if --pickle option is used. All '
              'directories must exists along the path.')
    )
    parser.add_argument(
        '-m', '--add-stmt-metadata',
        action='store_true',
        dest='add_stmt_metadata',
        help=('Optionally add special metadata to the evidence of all '
              'Statements that are produced, including the content IDs (the '
              'basenames of the files) and the readers used.')
    )
    parser.add_argument(
        '-p', '--pickle',
        action='store_true',
        dest='pickle',
        help='Select to use pickles instead of JSON for the dumps.'
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
    if path.isdir(args.input_file):
        file_list = [path.join(args.input_file, fname)
                     for fname in listdir(args.input_file)]
    elif path.isfile(args.input_file):
        with open(args.input_file, 'r') as f:
            file_list = [line.strip() for line in f.readlines()]
    else:
        raise ValueError("File or directory %s does not exist."
                         % args.input_file)

    logger.info("Found %d files." % len(file_list))
    for ftype in ['nxml', 'txt']:
        logger.debug('%d are %s' % (
            len([f for f in file_list if f.endswith(ftype)]), ftype
        ))

    # Select only a sample of the lines, if sample is chosen.
    if args.n_samp is not None:
        file_list = random.sample(file_list, args.n_samp)

    # If a range is specified, only use that range.
    if args.range_str is not None:
        start_idx, end_idx = [int(n) for n in args.range_str.split(':')]
        file_list = file_list[start_idx:end_idx]

    # Create a single base directory
    base_dir = get_dir('run_%s' % ('_and_'.join(args.readers)))

    # Set the verbosity. The quiet argument overrides the verbose argument.
    verbose = args.verbose and not args.quiet

    # Get the readers objects.
    readers = [reader_class(base_dir=base_dir, n_proc=args.n_proc)
               for reader_class in get_reader_classes()
               if reader_class.name.lower() in args.readers]

    # Read the files.
    outputs = read_files(file_list, readers, verbose=verbose)

    # Dump the outputs
    reading_out_path = path.join(args.output_path, 'readings')
    if args.pickle:
        reading_out_path += '.pkl'
        with open(reading_out_path, 'wb') as f:
            pickle.dump(outputs, f)
    else:
        reading_out_path += '.json'
        dump_readings(outputs, reading_out_path)
    print("Reading outputs stored in %s." % reading_out_path)

    # Generate and dump the statements.
    stmts_dump_path = path.join(args.output_path, 'stmts')
    stmt_gen = (s for rd in outputs for s in rd.get_statements())
    if args.pickle:
        stmts_dump_path += ".pkl"
        with open(stmts_dump_path, 'wb') as f:
            pickle.dump(list(stmt_gen), f)
    else:
        stmt_jsons = [s.to_json() for s in stmt_gen]
        stmts_dump_path += '.json'
        with open(stmts_dump_path, 'w') as f:
            json.dump(stmt_jsons, f)
    print("Statements stored in %s." % stmts_dump_path)


if __name__ == '__main__':
    main()
