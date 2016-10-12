from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import pickle

if __name__ == '__main__':
    file_list = sys.argv[1:]

    all_stmts = {}

    for file in file_list:
        with open(file, 'rb') as f:
            stmts = pickle.load(f)
        all_stmts.update(stmts)

    with open('reach_stmts.pkl', 'wb') as f:
        pickle.dump(all_stmts, f, protocol=2)

