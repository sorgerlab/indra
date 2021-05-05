import pickle
from collections import Counter
import numpy as np
from indra.belief import BeliefScorer


class SklearnScorer(BeliefScorer):
    """Uses a pre-trained Sklearn classifier to predict belief scores.

    Parameters
    ----------
    model_wrap : <indra.belief.sklearn.wrapper.SklearnBase>
        An instance of a Sklearn wrapper around an sklearn classifier.

    source_list : list of str
        Source API types in the same order that the model was trained on.
    prior_probs : dict
        Dictionary mapping source APIs (usually curated databases) that the
        model was not trained on to error parameters for the source. These
        error parameters are treated as independent.
    """
    def __init__(self, model_wrap):
        self.model_wrap = model_wrap

    """
    def score_statement(self, st, extra_evidence=None):
        belief_arr = self.score_statements[st]
        return belief_arr[0]
    """

    def check_prior_probs(self, statements):
        # TODO: Implement this
        return

    """
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
    """

    def score_statements(self, statements, extra_evidence=None):
        """Run predictions."""
        belief_arr = self.model_wrap.predict_proba(
                                statements, extra_evidence)[:, 1]
        return belief_arr
