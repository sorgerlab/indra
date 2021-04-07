import logging
from copy import deepcopy
from . import ModelChecker, NodesContainer
from indra.statements import *
from indra.ontology.bio import bio_ontology
from .model_checker import signed_edges_to_signed_nodes

logger = logging.getLogger(__name__)


class PybelModelChecker(ModelChecker):
    """Check a PyBEL model against a set of INDRA statements.

    Parameters
    ----------
    model : pybel.BELGraph
        A Pybel model to check.
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

    def get_graph(self, include_variants=False, symmetric_variant_links=False,
                  include_components=True, symmetric_component_links=True,
                  edge_filter_func=None):
        """Convert a PyBELGraph to a graph with signed nodes."""
        # This import is done here rather than at the top level to avoid
        # making pybel an implicit dependency of the model checker
        from indra.assemblers.pybel.assembler import belgraph_to_signed_graph
        if self.graph:
            return self.graph
        signed_edges = belgraph_to_signed_graph(
            self.model,
            include_variants=include_variants,
            symmetric_variant_links=symmetric_variant_links,
            include_components=include_components,
            symmetric_component_links=symmetric_component_links,
            propagate_annotations=True)
        self.graph = signed_edges_to_signed_nodes(
            signed_edges, copy_edge_data={'belief'})
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
        if obj is None:
            # Cannot check modifications for statements without object
            if isinstance(stmt, Modification):
                return (None, None, 'STATEMENT_TYPE_NOT_HANDLED')
            subj_nodes = self.get_nodes(subj, self.graph, 0)
            obj_nodes = self.get_nodes(obj, self.graph, 0)
        else:
            # Get the polarity for the statement
            if isinstance(stmt, Modification):
                target_polarity = 1 if isinstance(stmt, RemoveModification) \
                    else 0
                obj_agent = deepcopy(obj)
                obj_agent.mods.append(stmt._get_mod_condition())
                obj = obj_agent
            elif isinstance(stmt, RegulateActivity):
                target_polarity = 0 if stmt.is_activation else 1
                obj_agent = deepcopy(obj)
                obj_agent.activity = stmt._get_activity_condition()
                obj_agent.activity.is_active = True
                obj = obj_agent
            elif isinstance(stmt, RegulateAmount):
                target_polarity = 1 if isinstance(stmt, DecreaseAmount) else 0
            elif isinstance(stmt, Influence):
                target_polarity = 1 if stmt.overall_polarity() == -1 else 0

            subj_nodes = self.get_nodes(subj, self.graph, 0)
            obj_nodes = self.get_nodes(obj, self.graph, target_polarity)
            # Statement has object but it's not in the graph
            if not obj_nodes.all_nodes:
                return (None, None, 'OBJECT_NOT_FOUND')
            if not subj_nodes.all_nodes:
                return (None, None, 'SUBJECT_NOT_FOUND')
        return (subj_nodes, obj_nodes, None)

    def get_nodes(self, agent, graph, target_polarity):
        """Get all nodes corresponding to a given agent."""
        # This import is done here rather than at the top level to avoid
        # making pybel an implicit dependency of the model checker
        from indra.assemblers.pybel.assembler import _get_agent_node
        nc = NodesContainer(agent)
        if agent is None:
            nc.all_nodes = None
            return nc
        # First get exact match
        agent_node = _get_agent_node(agent)[0]
        if agent_node:
            node = (agent_node, target_polarity)
            if node in graph.nodes:
                nc.main_nodes.append(node)
        # Try get refined versions
        for n, ag in self.nodes_to_agents.items():
            if ag is not None and not ag.matches(agent) and ag.refinement_of(
                    agent, bio_ontology):
                node = (n, target_polarity)
                if node in graph.nodes:
                    nc.ref_nodes.append(node)
        nc.get_all_nodes()
        return nc

    def get_nodes_to_agents(self):
        """Return a dictionary mapping PyBEL nodes to INDRA agents."""
        if self.nodes_to_agents:
            return self.nodes_to_agents

        # This import is done here rather than at the top level to avoid
        # making pybel an implicit dependency of the model checker
        from indra.sources.bel.processor import get_agent
        self.nodes_to_agents = {
            node: get_agent(node) for node in self.model.nodes}
        return self.nodes_to_agents
