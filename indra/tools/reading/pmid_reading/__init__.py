from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import sys
import random
import shutil
import pickle
import logging
import argparse
from datetime import datetime
from platform import system
from .read_sparser import run_sparser
from .read_reach import run_reach, upload_process_reach_files

READER_DICT = {'reach': run_reach, 'sparser': run_sparser}


logger = logging.getLogger('runreader')


def get_proc_num():
    if system() == 'Linux':
        with open('/proc/cpuinfo', 'r') as f:
            ret = len([
                line for line in f.readlines() if line.startswith('processor')
                ])
    else:
        ret = None
    return ret


parser = argparse.ArgumentParser(
    description=('Apply NLP readers to the content available for a list of '
                 'pmids.')
    )
parser.add_argument(
    '-r', '--reader',
    choices=['reach', 'sparser', 'all'],
    default=['all'],
    dest='readers',
    nargs='+',
    help='Choose which reader(s) to use.'
    )
parser.add_argument(
    '-u', '--upload_json',
    dest='upload_json',
    action='store_true',
    help=('Option to simply upload previously read json files. Overrides -r '
          'option, so no reading will be done.')
    )
parser.add_argument(
    '--force_fulltext',
    dest='force_fulltext',
    action='store_true',
    help='Option to force reading of the full text.'
    )
parser.add_argument(
    '--force_read',
    dest='force_read',
    action='store_true',
    help='Option to force the reader to reread everything.'
    )
parser.add_argument(
    '-n', '--num_cores',
    dest='num_cores',
    default=1,
    type=int,
    help='Select the number of cores you want to use.'
    )
parser.add_argument(
    '-v', '--verbose',
    dest='verbose',
    action='store_true',
    help='Show more output to screen.'
    )
parser.add_argument(
    '-m', '--messy',
    dest='cleanup',
    action='store_false',
    help='Choose to not clean up after run.'
    )
parser.add_argument(
    '-s', '--start_index',
    dest='start_index',
    type=int,
    help='Select the first pmid in the list to start reading.',
    default=0
    )
parser.add_argument(
    '-e', '--end_index',
    dest='end_index',
    type=int,
    help='Select the last pmid in the list to read.',
    default=None
    )
parser.add_argument(
    '--shuffle',
    dest='shuffle',
    action='store_true',
    help=('Select a random sample of the pmids provided. -s/--start_index '
          'will be ingored, and -e/--end_index will set the number of '
          'samples to take.')
    )
parser.add_argument(
    '-o', '--outdir',
    dest='out_dir',
    default=None,
    help=('The output directory where stuff is written. This is only a '
          'temporary directory when reading. By default this will be the'
          '"<basename>_out".')
    )
parser.add_argument(
    dest='basename',
    help='The name of this job.'
    )
parser.add_argument(
    dest='pmid_list_file',
    help=('Path to a file containing a list of line separated pmids for the '
          'articles to be read.')
    )




if __name__ == '__main__':
    args = parser.parse_args()


def main(args):
    now = datetime.now()    # Set some variables
    if args.upload_json:
        args.readers = 'none'
    out_dir = args.out_dir
    if out_dir is None:
        out_dir = args.basename + '_out'
    made_outdir = False
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        made_outdir = True
    ret = None

    available_cores = get_proc_num()
    if args.num_cores >= available_cores:
        msg = ("You requested %d cores, but only %d available.\n" %
               (args.num_cores, available_cores) +
               "Are you sure you want to proceed? [y/N] > ")
        if sys.version_info.major > 2:
            resp = input(msg)
        else:
            resp = raw_input(msg)
        if resp.lower() not in ['y', 'yes']:
            logger.info("Aborting...")
            return
    try:
        # Option -u: just upload previously read JSON files
        if args.upload_json:
            with open(args.pmid_list_file, 'rb') as f:
                text_sources = pickle.load(f)
            stmts = upload_process_reach_files(
                args.out_dir,
                text_sources,
                args.num_cores
                )
            pickle_file = '%s_stmts.pkl' % args.basename
            with open(pickle_file, 'wb') as f:
                pickle.dump(stmts, f, protocol=2)
            sys.exit()

        # Option -r <reader>: actually read the content.

        # Load the list of PMIDs from the given file
        with open(args.pmid_list_file) as f:
            pmid_list = [line.strip('\n') for line in f.readlines()]

        if args.end_index is None:
            args.end_index = len(pmid_list)

        if args.shuffle:
            pmid_list = random.sample(pmid_list, args.end_index)

        # Do the reading
        readers = []
        if 'all' in args.readers:
            readers = list(READER_DICT.keys())
        else:
            readers = args.readers[:]

        stmts = {}
        for reader in readers:
            run_reader = READER_DICT[reader]
            some_stmts, _ = run_reader(
                pmid_list,
                out_dir,
                args.num_cores,
                args.start_index,
                args.end_index,
                args.force_read,
                args.force_fulltext,
                cleanup=args.cleanup,
                verbose=args.verbose
                )
            stmts[reader] = some_stmts

        N_tot = sum([
            len(stmts[reader][pmid]) for reader in readers
            for pmid in stmts[reader].keys()
            ])
        logger.info('Collected a total of %s statements.' % N_tot)
        for reader in readers:
            N = sum([len(l) for l in stmts[reader].values()])
            logger.info(
                '%s accumulated %d statements.' % (reader.capitalize(), N)
                )

        # Pickle the statements
        if args.end_index is None:
            args.end_index = 'end'
        pickle_file = '%s_stmts_%s-%s.pkl' % (
            args.basename,
            args.start_index,
            args.end_index
            )
        with open(pickle_file, 'wb') as f:
            pickle.dump(stmts, f, protocol=2)
        ret = pickle_file
    finally:
        time_taken = datetime.now() - now
        print('This run took', time_taken)
        time_file = os.path.join(os.path.dirname(__file__), 'time_data.txt')
        with open(time_file, 'a') as f:
            f.write('Started run at %s with args %s lasting %s.\n' %
                    (now, str(args), time_taken))
        if made_outdir and args.cleanup:
            shutil.rmtree(out_dir)
    return ret


if __name__ == '__main__':
    main(args)
