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
                      'belief_sklearn_test_stmts.pkl')

test_df_path = join(dirname(abspath(__file__)),
                    'belief_sklearn_test_df.pkl')

with open(test_stmt_path, 'rb') as f:
    test_stmts, y_arr_stmts = pickle.load(f)

with open(test_df_path, 'rb') as f:
    test_df, y_arr_df = pickle.load(f)


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


# TODO: Made this so it's not a ValueError, this may change back in the future
#@raises(ValueError)
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


def test_fit_stmts():
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'signor']
    cw = CountsModel(lr, source_list)
    cw.fit(test_stmts, y_arr_stmts)
    # Once the model is fit, the coef_ attribute should be defined
    assert 'coef_' in cw.model.__dict__


def test_fit_stmts_predict_stmts():
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'signor']
    cw = CountsModel(lr, source_list)
    cw.fit(test_stmts, y_arr_stmts)
    probs = cw.predict_proba(test_stmts)
    assert probs.shape == (len(test_stmts), 2), \
        'prediction results should have dimension (# stmts, # classes)'
    log_probs = cw.predict_log_proba(test_stmts)
    assert log_probs.shape == (len(test_stmts), 2), \
        'prediction results should have dimension (# stmts, # classes)'
    preds = cw.predict(test_stmts)
    assert preds.shape == (len(test_stmts),), \
        'prediction results should have dimension (# stmts)'


@raises(ValueError)
def test_check_df_cols_err():
    """Drop a required column and make sure we get a ValueError."""
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'signor']
    cw = CountsModel(lr, source_list)
    cw.df_to_matrix(test_df.drop('agB_ns', axis=1))


def test_check_df_cols_noerr():
    """Test dataframe should not raise ValueError."""
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'signor']
    cw = CountsModel(lr, source_list)
    cw.df_to_matrix(test_df)


def test_df_to_matrix():
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'signor']
    cw = CountsModel(lr, source_list)
    x_arr = cw.df_to_matrix(test_df)
    #import ipdb; ipdb.set_trace()
    assert isinstance(x_arr, np.ndarray), 'x_arr should be a numpy array'
    assert x_arr.shape == (len(test_df), len(source_list)), \
            'stmt matrix dimensions should match test stmts'
    assert x_arr.shape == (len(test_df), len(source_list))
    # Try again with statement type
    cw = CountsModel(lr, source_list, use_stmt_type=True)
    x_arr = cw.df_to_matrix(test_df)
    assert x_arr.shape == (len(test_df), len(source_list)+1), \
        'matrix should have a col for sources and stmt type'


def test_fit_df():
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'medscan', 'trips', 'rlimsp']
    cw = CountsModel(lr, source_list)
    cw.fit(test_df, y_arr_df)
    # Once the model is fit, the coef_ attribute should be defined
    assert 'coef_' in cw.model.__dict__


def test_fit_stmts_pred_df():
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'signor']
    cw = CountsModel(lr, source_list)
    # Train on statement data
    cw.fit(test_stmts, y_arr_stmts)
    # Predict on DF data
    probs = cw.predict_proba(test_df)
    assert probs.shape == (len(test_df), 2), \
        'prediction results should have dimension (# stmts, # classes)'
    log_probs = cw.predict_log_proba(test_df)
    assert log_probs.shape == (len(test_df), 2), \
        'prediction results should have dimension (# stmts, # classes)'
    preds = cw.predict(test_df)
    assert preds.shape == (len(test_df),), \
        'prediction results should have dimension (# stmts)'


def test_fit_df_pred_stmts():
    lr = LogisticRegression()
    source_list = ['reach', 'sparser', 'signor']
    cw = CountsModel(lr, source_list)
    # Train on statement data
    cw.fit(test_df, y_arr_df)
    # Predict on DF data
    probs = cw.predict_proba(test_stmts)
    assert probs.shape == (len(test_stmts), 2), \
        'prediction results should have dimension (# stmts, # classes)'
    log_probs = cw.predict_log_proba(test_stmts)
    assert log_probs.shape == (len(test_stmts), 2), \
        'prediction results should have dimension (# stmts, # classes)'
    preds = cw.predict(test_stmts)
    assert preds.shape == (len(test_stmts),), \
        'prediction results should have dimension (# stmts)'


if __name__ == '__main__':
    test_fit_stmts_pred_df()
