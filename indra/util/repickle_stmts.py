import sys
import pickle

filename = sys.argv[1]

with open(filename, 'rb') as f:
    stmts = pickle.load(f)

with open('%s.new' % filename, 'wb') as f:
    pickle.dump(stmts, f)

