import logging
import networkx as nx
from . import ModelChecker, NodesContainer
from indra.statements import *
from indra.ontology.bio import bio_ontology


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
        self.graph = nx.DiGraph()
        nodes = []
        for node, node_data in self.model.nodes(data=True):
            nodes.append(((node, 0), node_data))
        self.graph.add_nodes_from(nodes)
        for (u, v, data) in self.model.edges(data=True):
            self.graph.add_edge((u, 0), (v, 0), belief=data['belief'])
        self.get_nodes_to_agents()
        return self.graph

    def process_statement(self, stmt):
        # Check if this is one of the statement types that we can check
        if not isinstance(stmt, (Modification, RegulateAmount,
                                 RegulateActivity, Influence)):
            logger.info('Statement type %s not handled' %
                        stmt.__class__.__name__)
            return (None, None, 'STATEMENT_TYPE_NOT_HANDLED')
        subj, obj = stmt.agent_list()
        subj_nodes = self.get_nodes(subj, self.graph)
        obj_nodes = self.get_nodes(obj, self.graph)
        # Statement has object but it's not in the graph
        if obj and not obj_nodes.all_nodes:
            return (None, None, 'OBJECT_NOT_FOUND')
        if subj and not subj_nodes.all_nodes:
            return (None, None, 'SUBJECT_NOT_FOUND')
        return (subj_nodes, obj_nodes, None)

    def _sample_paths(self, input_set, obj_name, target_polarity,
                      max_paths=1, max_path_length=5):
        # TODO implement sampling
        pass

    def get_nodes(self, agent, graph):
        """Get all nodes corresponding to a given agent."""
        nc = NodesContainer(agent)
        if agent is None:
            nc.all_nodes = None
            return nc
        node = (agent.name, 0)
        if node in graph.nodes:
            nc.main_nodes.append(node)
        for n, ag in self.nodes_to_agents.items():
            if ag is not None and not ag.matches(agent) and ag.refinement_of(
                    agent, bio_ontology):
                node = (n, 0)
                if node in graph.nodes:
                    nc.ref_nodes.append(node)
        nc.get_all_nodes()
        return nc

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
