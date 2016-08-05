import pickle

file_list = [
	'reach_stmts_0_1500.pkl',
	'reach_stmts_1500_2500.pkl',
	'reach_stmts_2500_3500.pkl',
	'reach_stmts_3500_4500.pkl',
	'reach_stmts_4500_5500.pkl',
	'reach_stmts_5500_6671.pkl',
]

all_stmts = {}

for file in file_list:
    with open(file) as f:
        stmts = pickle.load(f)
    all_stmts.update(stmts)

with open('reach_stmts.pkl', 'w') as f:
    pickle.dump(all_stmts, f)

