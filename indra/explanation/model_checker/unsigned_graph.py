import logging
import networkx as nx
from . import ModelChecker
from indra.statements import *


logger = logging.getLogger(__name__)


class UnsignedGraphModelChecker(ModelChecker):
    """Check an unsigned DiGraph against a set of INDRA statements.

    Parameters
    ----------
    model : networkx.DiGraph
        Unsigned DiGraph to check.
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
        self.graph = nx.DiGraph()
        nodes = []
        for node, node_data in self.model.nodes(data=True):
            nodes.append(((node, 0), node_data))
        self.graph.add_nodes_from(nodes)
        for (u, v, data) in self.model.edges(data=True):
            self.graph.add_edge((u, 0), (v, 0), belief=data['belief'])
        return self.graph

    def process_statement(self, stmt):
        # Check if this is one of the statement types that we can check
        if not isinstance(stmt, (Modification, RegulateAmount,
                                 RegulateActivity, Influence)):
            logger.info('Statement type %s not handled' %
                        stmt.__class__.__name__)
            return (None, None, 'STATEMENT_TYPE_NOT_HANDLED')
        subj, obj = stmt.agent_list()
        if subj is None or (subj.name, 0) not in self.graph.nodes:
            return (None, None, 'SUBJECT_NOT_FOUND')
        if obj is None or (obj.name, 0) not in self.graph.nodes:
            return (None, None, 'OBJECT_NOT_FOUND')
        return ([(subj.name, 0)], [(obj.name, 0)], None)

    def process_subject(self, subj):
        return [subj], None

    def _sample_paths(self, input_set, obj_name, target_polarity,
                      max_paths=1, max_path_length=5):
        # TODO implement sampling
        pass
