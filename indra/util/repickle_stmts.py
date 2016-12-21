from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

if __name__ == '__main__':
    import sys
    import pickle

    filename = sys.argv[1]

    with open(filename, 'rb') as f:
        stmts = pickle.load(f)

    with open('%s.new' % filename, 'wb') as f:
        pickle.dump(stmts, f, protocol=2)
