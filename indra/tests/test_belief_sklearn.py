import random
import pickle
import numpy as np
from collections import defaultdict
from os.path import join, abspath, dirname
from nose.tools import raises
from sklearn.linear_model import LogisticRegression
from indra.sources import signor
from indra.belief.sklearn.wrapper import CountsModel

test_stmt_path = join(dirname(abspath(__file__)),
                      'belief_sklearn_test_data.pkl')
with open(test_stmt_path, 'rb') as f:
    test_stmts, y_arr = pickle.load(f)

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
    source_list = ['reach', 'sparser']
    cw = CountsModel(lr, source_list)

@raises(ValueError)
def test_missing_source():
    """Check that all source_apis in training data are in source list."""
    lr = LogisticRegression()
    source_list = ['reach', 'sparser']
    cw = CountsModel(lr, source_list)
    # Should error because test stmts are from signor and signor
    # is not in list
    cw.stmts_to_matrix(test_stmts)

def test_stmts_to_matrix():
    """Check that all source_apis in training data are in source list."""
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'signor']
    cw = CountsModel(lr, source_list)
    x_arr = cw.stmts_to_matrix(test_stmts)
    assert isinstance(x_arr, np.ndarray), 'x_arr should be a numpy array'
    assert x_arr.shape == (len(test_stmts), len(source_list)), \
            'stmt matrix dimensions should match test stmts'
    assert set(x_arr.sum(axis=0)) == set([0, 0, len(test_stmts)]), \
           'Signor col should be 1 in every row, other cols 0.'
    # Try again with statement type
    cw = CountsModel(lr, source_list, use_stmt_type=True)
    x_arr = cw.stmts_to_matrix(test_stmts)
    assert x_arr.shape == (len(test_stmts), len(source_list)+1), \
        'matrix should have a col for sources and stmt type'

