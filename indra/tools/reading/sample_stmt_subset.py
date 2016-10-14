from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

if __name__ == '__main__':
    import sys
    import numpy as np
    import pickle

    if len(sys.argv) != 3:
        print("Usage: %s stmt_file num_sample_papers" % sys.argv[0])
        sys.exit()

    input_file = sys.argv[1]
    sample_size = int(sys.argv[2])
    print("Opening input file...")
    with open(input_file, 'rb') as f:
        stmts_by_paper = pickle.load(f)

    pmids = stmts_by_paper.keys()
    if sample_size >= len(pmids):
        print("Sample size exceeds the total number of papers.")
        print("No need to sample.")
        sys.exit()

    sample_ids = np.random.choice(pmids, sample_size, replace=False)
    print("Building subset...")
    subset = {paper_id: stmts_by_paper[paper_id] for paper_id in sample_ids}

    print("Pickling paper subset...")
    with open('%s_subset.pkl' % sys.argv[1], 'wb') as f:
        pickle.dump(subset, f, protocol=2)

