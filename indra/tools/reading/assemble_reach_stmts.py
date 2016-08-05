import sys
import pickle

file_list = sys.argv[1:]

all_stmts = {}

for file in file_list:
    with open(file) as f:
        stmts = pickle.load(f)
    all_stmts.update(stmts)

with open('reach_stmts.pkl', 'w') as f:
    pickle.dump(all_stmts, f)

