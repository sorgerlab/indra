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

    def score_statements(self, statements, extra_evidence=None):
        """Run predictions."""
        belief_arr = self.model_wrap.predict_proba(
                                statements, extra_evidence)[:, 1]
        return belief_arr
