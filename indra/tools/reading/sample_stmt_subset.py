import sys
import numpy as np
import pickle

if __name__ == '__main__':

    if len(sys.argv) != 3:
        print "Usage: %s stmt_file num_sample_papers" % sys.argv[0]
        sys.exit()

    input_file = sys.argv[1]
    sample_size = int(sys.argv[2])
    print "Opening input file..."
    with open(input_file) as f:
        stmts_by_paper = pickle.load(f)

    pmids = stmts_by_paper.keys()
    sample_ids = np.random.choice(pmids, sample_size, replace=False)
    print "Building subset..."
    subset = {paper_id: stmts_by_paper[paper_id] for paper_id in sample_ids}

    print "Pickling paper subset..."
    with open('%s_subset.pkl' % sys.argv[1], 'w') as f:
        pickle.dump(subset, f)

