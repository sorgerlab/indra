import logging
from . import ModelChecker
from indra.statements import *
from ..pathfinding_util import signed_edges_to_signed_nodes

logger = logging.getLogger(__name__)


class SignedGraphModelChecker(ModelChecker):
    """Check an signed MultiDiGraph against a set of INDRA statements.

    Parameters
    ----------
    model : networkx.MultiDiGraph
        Signed MultiDiGraph to check.
    statements : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to check the model against.
    do_sampling : bool
        Whether to use breadth-first search or weighted sampling to
        generate paths. Default is False (breadth-first search).
    seed : int
        Random seed for sampling (optional, default is None).
    """
    def __init__(self, model, statements=None, do_sampling=False, seed=None):
        super().__init__(model, statements, do_sampling, seed)

    def get_graph(self):
        if self.graph:
            return self.graph
        self.graph = signed_edges_to_signed_nodes(self.model)
        return self.graph

    def process_statement(self, stmt):
        # Check if this is one of the statement types that we can check
        if not isinstance(stmt, (Activation, Inhibition,
                                 IncreaseAmount, DecreaseAmount, Influence)):
            logger.info('Statement type %s not handled' %
                        stmt.__class__.__name__)
            return (None, None, 'STATEMENT_TYPE_NOT_HANDLED')
        # Get the polarity for the statement
        if isinstance(stmt, RegulateActivity):
            target_polarity = 0 if stmt.is_activation else 1
        elif isinstance(stmt, RegulateAmount):
            target_polarity = 1 if isinstance(stmt, DecreaseAmount) else 0
        elif isinstance(stmt, Influence):
            target_polarity = 1 if stmt.overall_polarity() == -1 else 0
        subj, obj = stmt.agent_list()
        if subj is None or (subj.name, 0) not in self.graph.nodes:
            return (None, None, 'SUBJECT_NOT_FOUND')
        if obj is None or (obj.name, target_polarity) not in self.graph.nodes:
            return (None, None, 'OBJECT_NOT_FOUND')
        return ([(subj.name, 0)], [(obj.name, target_polarity)], None)

    def process_subject(self, subj):
        return [subj], None
