import logging
from . import ModelChecker
from indra.statements import *
from indra.ontology.bio import bio_ontology
from .model_checker import signed_edges_to_signed_nodes

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
    nodes_to_agents : dict
        A dictionary mapping nodes of intermediate signed edges graph to INDRA
        agents.

    Attributes
    ----------
    graph : nx.Digraph
        A DiGraph with signed nodes to find paths in.
    """
    def __init__(self, model, statements=None, do_sampling=False, seed=None,
                 nodes_to_agents=None):
        super().__init__(model, statements, do_sampling, seed, nodes_to_agents)

    def get_graph(self):
        if self.graph:
            return self.graph
        self.graph = signed_edges_to_signed_nodes(
            self.model, copy_edge_data={'belief'})
        self.get_nodes_to_agents()
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
        if obj is None:
            obj_nodes = [None]
        else:
            obj_nodes = self.get_nodes(obj, self.graph, target_polarity)
            # Statement has object but it's not in the graph
            if not obj_nodes:
                return (None, None, 'OBJECT_NOT_FOUND')
        return ([subj], obj_nodes, None)

    def process_subject(self, subj):
        # We will not get here if subject is None
        subj_nodes = self.get_nodes(subj, self.graph, 0)
        # Statement has subject but it's not in the graph
        if not subj_nodes:
            return (None, 'SUBJECT_NOT_FOUND')
        return subj_nodes, None

    def get_nodes(self, agent, graph, target_polarity):
        """Get all nodes corresponding to a given agent."""
        nodes = set()
        node = (agent.name, target_polarity)
        if node in graph.nodes:
            nodes.add(node)
        for n, ag in self.nodes_to_agents.items():
            if ag is not None and ag.refinement_of(agent, bio_ontology):
                node = (n, target_polarity)
                if node in graph.nodes:
                    nodes.add(node)
        return nodes

    def get_nodes_to_agents(self):
        """Return a dictionary mapping IndraNet nodes to INDRA agents."""
        if self.nodes_to_agents:
            return self.nodes_to_agents

        # NOTE: this way of retrieving agents might miss some important
        # agent properties. The recommended way is to provide this mapping
        # externally.
        graph = self.get_graph()
        for node, data in graph.nodes(data=True):
            ag = Agent(node[0], db_refs={data.get('ns'): data.get('id')})
            self.nodes_to_agents[node[0]] = ag
        return self.nodes_to_agents
