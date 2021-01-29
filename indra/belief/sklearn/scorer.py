import pickle
from collections import Counter
import numpy as np
from indra.belief import BeliefScorer


class SklearnScorer(BeliefScorer):
    """Uses a pre-trained Sklearn classifier to predict belief scores.

    Parameters
    ----------
    model_pkl : str
        Pickle file containing the sklearn classifier.
    source_list : list of str
        Source API types in the same order that the model was trained on.
    prior_probs : dict
        Dictionary mapping source APIs (usually curated databases) that the
        model was not trained on to error parameters for the source. These
        error parameters are treated as independent.
    """
    def __init__(self, model):
        self.model = model
        #self.source_list = source_list
        #self.correct_ix = correct_ix
        # Load the model from file
        #with open(self.model_pkl, 'rb') as f:
        #    self.model = pickle.load(f)

    def score_statement(self, st, extra_evidence=None):
        # Combine the evidences
        if extra_evidence is None:
            extra_evidence = []
        evidences = st.evidence + extra_evidence
        # Just in case there are no evidences...
        if not evidences:
            return 0
        # Collect all unique sources
        sources = [ev.source_api for ev in evidences]
        src_ctr = Counter(sources)
        # Numpy array of counts in correct order
        counts = np.array([[src_ctr.get(src, 0) for src in self.source_list]])
        # Predict!
        belief = self.model.predict_proba(counts)[0][self.correct_ix]
        return belief

    def score_statements(self, statements):
        """Run predictions."""
        # Build training data matrix
        num_cols = len(self.source_list)
        num_rows = len(statements)
        arr = np.zeros((num_rows, num_cols))
        for stmt_ix, stmt in enumerate(statements):
            sources = [ev.source_api for ev in stmt.evidence]
            src_ctr = Counter(sources)
            for src_ix, src in enumerate(self.source_list):
                arr[stmt_ix, src_ix] = src_ctr.get(src, 0)
        scores = self.model.predict_proba(arr)
        for stmt_ix, stmt in enumerate(statements):
            stmt.belief = scores[stmt_ix, self.correct_ix]
