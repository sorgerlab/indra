import random
import pickle
from collections import defaultdict
from os.path import join, abspath, dirname
from sklearn.linear_model import LogisticRegression
from indra.sources import signor
from indra.belief.sklearn.wrapper import CountsModel

test_stmt_path = join(dirname(abspath(__file__)),
                      'belief_sklearn_test_data.pkl')
#with open(test_stmt_path, 'rb') as f:
#    test_stmts, y_arr = pickle.load(f)

# A set of statements derived from Signor used for testing purposes.
def _dump_test_data(filename, num_per_type=10):
    """Get corpus of statements for testing that has a range of stmt types."""
    sp = signor.process_from_web()
    # Group statements by type
    stmts_by_type = defaultdict(list)
    for stmt in sp.statements:
        stmts_by_type[stmt.__class__].append(stmt)
    # Sample statements of each type (without replacement)
    stmt_sample = []
    for stmt_type, stmt_list in stmts_by_type.items():
        if len(stmt_list) <= num_per_type:
            stmt_sample.extend(stmt_list)
        else:
            stmt_sample.extend(random.sample(stmt_list, num_per_type))
    # Make a random binary class vector for the stmt list
    y_arr = [random.choice((0, 1)) for s in stmt_sample]
    with open(test_stmt_path, 'wb') as f:
        pickle.dump((stmt_sample, y_arr), f)
    return stmt_sample

def test_counts_wrapper():
    """Instantiate counts wrapper and make stmt matrix"""
    lr = LogisticRegression()
    cw = CountsModel(lr)
