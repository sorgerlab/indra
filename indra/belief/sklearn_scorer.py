import pickle
from collections import Counter
import numpy as np
from . import BeliefScorer

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
    def __init__(self, model_pkl, source_list, correct_ix=1):
        self.model_pkl = model_pkl
        self.source_list = source_list
        self.correct_ix = correct_ix
        # Load the model from file
        with open(self.model_pkl, 'rb') as f:
            self.model = pickle.load(f)

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

