from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from argparse import ArgumentParser

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Assemble many pickle files into one.'
        )
    parser.add_argument(
        '-r', '--readers',
        dest='readers',
        nargs='+',
        help='Choose which reader(s) to use.'
        )
    parser.add_argument(
        dest='file_list',
        nargs='+',
        help='A list of file paths.'
        )
    args = parser.parse_args()

    all_stmts = dict.fromkeys(args.readers)
    for k in all_stmts.keys():
        all_stmts[k] = {}

    for file in args.file_list:
        with open(file, 'rb') as f:
            stmts = pickle.load(f)
        for reader in args.readers:
            all_stmts[reader].update(stmts[reader])

    with open('reading_stmts.pkl', 'wb') as f:
        pickle.dump(all_stmts, f, protocol=2)
