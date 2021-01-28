from collections import Counter
import numpy as np
from sklearn.preprocessing import OrdinalEncoder
from indra.statements import get_all_descendants, Statement

class SklearnBase(object):
    """Base class to wrap an Sklearn model with statement preprocessing.

    Parameters
    ----------
    model : sklearn or similar model
        Any instance of a classifier object supporting the methods `fit`
        and `predict_proba`.
    """
    def __init__(self, model):
        self.model = model

    def stmts_to_matrix(self, stmts, y_arr, *args, **kwargs):
        return NotImplementedError('Need to implement the stmts_to_matrix '
                                   'method')

    def fit(self, stmts, y_arr, *args, **kwargs):
        # Preprocess the statements into features
        # TODO
        pass

    def predict_proba(self, stmts):
        stmts_arr = self.stmts_to_matrix(stmts)
        return self.model.predict_proba(stmts_arr)

class CountsModel(SklearnBase):
    """Predictor based on source evidence counts and other stmt properties.

    Parameters
    ----------
    source_list : list of str
        List of strings denoting the evidence sources (evidence.source_api
        values) used for prediction.
    use_stmt_type : bool
        Whether to include statement type as a feature.
    use_num_members : bool
        Whether have a feature denoting the number of members of the statement.
        Primarily for stratifying belief predictions about Complex statements
        with more than two members.
    """
    def __init__(self, model, source_list, use_stmt_type=False,
                 use_num_members=False):
        # Call superclass constructor to store the model
        super(CountsModel, self).__init__(model)
        self.use_stmt_type = use_stmt_type
        self.use_num_members = use_num_members
        self.source_list = source_list

        # Build dictionary mapping INDRA Statement types to integers
        if use_stmt_type:
            all_stmt_types = get_all_descendants(Statement)
            self.stmt_type_map = {t: ix for ix, t in enumerate(all_stmt_types)}


    def stmts_to_matrix(self, stmts, y_arr):
        # Check dimensions of stmts (x) and y_arr
        if len(stmts) != len(y_arr):
            raise ValueError("Number of stmts must match length of y_arr.")
        # Initialize a Numpy Array to store statement features and class
        # variable
        # Option 1: include stmt_type (encoded) and score
        #df = kge_join.drop(columns=['stmt_hash', 'stmt_num', 'agA_name',
        # 'agB_name'])

        # Add categorical features
        cat_features = []
        for stmt in stmts:
            feature_row = []
            if self.use_stmt_type:
                feature_row.append(self.stmt_type_map[type(stmt)])
            # Only add a feature row if we're using some of the features.
            if feature_row:
                cat_features.append(feature_row)

        # Get source count features
        num_cols = len(self.source_list)
        num_rows = len(stmts)
        x_arr = np.zeros((num_rows, num_cols))
        for stmt_ix, stmt in enumerate(stmts):
            sources = [ev.source_api for ev in stmt.evidence]
            src_ctr = Counter(sources)
            for src_ix, src in enumerate(self.source_list):
                arr[stmt_ix, src_ix] = src_ctr.get(src, 0)
        return arr

