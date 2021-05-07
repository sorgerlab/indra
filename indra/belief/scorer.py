import pickle
from collections import Counter
from typing import Sequence, Optional, List
import numpy as np
from indra.belief import BeliefScorer
from indra.statements import Evidence, Statement
from .wrapper import SklearnBase


class SklearnScorer(BeliefScorer):
    """Uses a pre-trained Sklearn classifier to predict belief scores.

    Parameters
    ----------
    model_wrap :
        An instance of a Sklearn wrapper around an sklearn classifier.
    """
    def __init__(
        self,
        model_wrap: SklearnBase
    ):
        self.model_wrap = model_wrap

    def check_prior_probs(
        self,
        statements: Sequence[Statement],
    ) -> None:
        """Empty implementation for now."""
        pass

    def score_statements(
        self,
        statements: Sequence[Statement],
        extra_evidence: Optional[List[List[Evidence]]] = None,
    ) -> Sequence[float]:
        """Run predictions."""
        belief_arr = self.model_wrap.predict_proba(
                                statements, extra_evidence)[:, 1]
        return belief_arr
