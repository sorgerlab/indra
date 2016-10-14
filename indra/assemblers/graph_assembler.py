from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import itertools
from indra.statements import *
import logging

logger = logging.getLogger('graph_assembler')
try:
    import pygraphviz
except ImportError:
    logger.error('Cannot use graph assembler because '
                 'pygraphviz could not be imported.')

default_graph_properties = {
    'directed': True,
    'fixedsize': True,
    'fontname': 'arial',
    'splines': 'spline',
    'rankdir': 'LR'
    }

default_node_properties = {
    'color': '#FBAF3F',
    'shape': 'Mrecord',
    'fontsize': 8
    }

default_edge_properties = {
    'arrowsize': 0.5
    }


class GraphAssembler():
    """The Graph assembler assembles INDRA Statements into a
    Graphviz node-edge graph.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be added to the assembler's list
        of Statements.
    graph_properties : Optional[dict[str: str]]
        A dictionary of graphviz graph properties overriding the default ones.
    node_properties : Optional[dict[str: str]]
        A dictionary of graphviz node properties overriding the default ones.
    edge_properties : Optional[dict[str: str]]
        A dictionary of graphviz edge properties overriding the default ones.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to be assembled.
    graph : pygraphviz.AGraph
        A pygraphviz graph that is assembled by this assembler.
    existing_nodes : list[tuple]
        The list of nodes (identified by node key tuples) that are
        alredy in the graph.
    existing_edges : list[tuple]
        The list of edges (identified by edge key tuples) that are
        already in the graph.
    graph_properties : dict[str: str]
        A dictionary of graphviz graph properties used for assembly.
    node_properties : dict[str: str]
        A dictionary of graphviz node properties used for assembly.
    edge_properties : dict[str: str]
        A dictionary of graphviz edge properties used for assembly.
        Note that most edge properties are determined based on the type of
        the edge by the assembler (e.g. color, arrowhead).
        These settings cannot be directly controlled through the API.
    """
    def __init__(self, stmts=None, graph_properties=None,
                 node_properties=None, edge_properties=None):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts
        self.graph_properties = default_graph_properties
        self.node_properties = default_node_properties
        self.edge_properties = default_edge_properties
        if graph_properties:
            for k, v in graph_properties.items():
                self.graph_properties[k] = v
        if node_properties:
            for k, v in node_properties.items():
                self.node_properties[k] = v
        if edge_properties:
            for k, v in edge_properties.items():
                self.edge_properties[k] = v
        self.graph = pygraphviz.AGraph(**self.graph_properties)
        self.existing_nodes = []
        self.existing_edges = []
        self._complex_nodes = []

    def add_statements(self, stmts):
        """Add a list of statements to be assembled.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of INDRA Statements to be appended to the assembler's list.
        """
        for stmt in stmts:
            self.statements.append(stmt)

    def make_model(self):
        """Assemble the graph from the assembler's list of INDRA Statements."""
        # Assemble in two stages.
        # First, create the nodes of the graph
        for stmt in self.statements:
            if isinstance(stmt, Phosphorylation):
                if stmt.enz is None:
                    continue
                self._add_node(stmt.enz)
                self._add_node(stmt.sub)
            elif isinstance(stmt, Dephosphorylation):
                if stmt.enz is None:
                    continue
                self._add_node(stmt.enz)
                self._add_node(stmt.sub)
            elif isinstance(stmt, Activation):
                self._add_node(stmt.subj)
                self._add_node(stmt.obj)
            elif isinstance(stmt, Complex):
                for m in stmt.members:
                    self._add_node(m)
        # Second, create the edges of the graph
        for stmt in self.statements:
            if isinstance(stmt, Phosphorylation):
                if stmt.enz is None:
                    continue
                self._add_phosphorylation(stmt.enz, stmt.sub)
            elif isinstance(stmt, Dephosphorylation):
                if stmt.enz is None:
                    continue
                self._add_dephosphorylation(stmt.enz, stmt.sub)
            elif isinstance(stmt, Activation):
                self._add_activation(stmt.subj, stmt.obj,
                                     stmt.is_activation)
            elif isinstance(stmt, Complex):
                self._add_complex(stmt.members)

    def get_string(self):
        """Return the assembled graph as a string.

        Returns
        -------
        graph_string : str
            The assembled graph as a string.
        """
        graph_string = self.graph.to_string()
        graph_string = graph_string.replace('\\N', '\\n')
        return graph_string

    def save_dot(self, file_name='graph.dot'):
        """Save the graph in a graphviz dot file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the graph dot string to.
        """
        s = self.get_string()
        with open(file_name, 'wt') as fh:
            fh.write(s)

    def save_pdf(self, file_name='graph.pdf', prog='dot'):
        """Draw the graph and save as an image or pdf file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the graph as. Default: graph.pdf
        prog : Optional[str]
            The graphviz program to use for graph layout. Default: dot
        """
        self.graph.draw(file_name, prog=prog)

    def _add_edge(self, source, target, **kwargs):
        """Add an edge to the graph."""
        # Start with default edge properties
        edge_properties = self.edge_properties
        # Overwrite ones that are given in function call explicitly
        for k, v in kwargs.items():
            edge_properties[k] = v
        self.graph.add_edge(source, target, **edge_properties)

    def _add_node(self, agent):
        """Add an Agent as a node to the graph."""
        if agent is None:
            return
        if not agent.bound_conditions:
            node_label = _get_node_label(agent)
        else:
            bound_agents = [bc.agent for bc in agent.bound_conditions if
                            bc.is_bound]
            if bound_agents:
                bound_names = [_get_node_label(a) for a in bound_agents]
                node_label = _get_node_label(agent) + '/' + \
                             '/'.join(bound_names)
                self._complex_nodes.append([agent] + bound_agents)
            else:
                node_label = _get_node_label(agent)
        node_key = _get_node_key(agent)
        if node_key in self.existing_nodes:
            return
        self.existing_nodes.append(node_key)
        self.graph.add_node(node_key,
                        label=node_label,
                        **self.node_properties)

    def _add_phosphorylation(self, enz, sub):
        """Assemble a Phosphorylation statement."""
        source = _get_node_key(enz)
        target = _get_node_key(sub)
        edge_key = (source, target, 'phosphorylation')
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        params = {'color': '#000000',
                  'arrowhead': 'normal',
                  'dir': 'forward'}
        self._add_edge(source, target, **params)

    def _add_dephosphorylation(self, enz, sub):
        """Assemble a Dephosphorylation statement."""
        source = _get_node_key(enz)
        target = _get_node_key(sub)
        edge_key = (source, target, 'dephosphorylation')
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        params = {'color': '#ff0000',
                  'arrowhead': 'normal',
                  'dir': 'forward'}
        self._add_edge(source, target, **params)

    def _add_activation(self, subj, obj, rel):
        """Assemble an Activation statment."""
        source = _get_node_key(subj)
        target = _get_node_key(obj)
        edge_key = (source, target, 'activation', rel)
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        color = '#000000' if rel else '#ff0000'
        arrowhead = 'vee' if rel else 'tee'
        params = {'color': color,
                  'arrowhead': arrowhead,
                  'dir': 'fowrard'}
        self._add_edge(source, target, **params)

    def _add_complex(self, members):
        """Assemble a Complex statement."""
        params = {'color': '#0000ff',
                  'arrowhead': 'dot',
                  'arrowtail': 'dot',
                  'dir': 'both'}
        for m1, m2 in itertools.combinations(members, 2):
            if self._has_complex_node(m1, m2):
                continue
            m1_key = _get_node_key(m1)
            m2_key = _get_node_key(m2)
            edge_key = (set([m1_key, m2_key]), 'complex')
            if edge_key in self.existing_edges:
                return
            self.existing_edges.append(edge_key)
            self._add_edge(m1_key, m2_key, **params)

    def _has_complex_node(self, m1, m2):
        for cplx in self._complex_nodes:
            names = [m.name for m in cplx]
            if m1.name in names and m2.name in names:
                return True
            else:
                return False

def _get_node_label(agent):
    # If the agent doesn't have grounding in a known
    # database, try to use the original text as a node name.
    # otherwise return the agent name.
    if not (agent.db_refs.get('UP') or
            agent.db_refs.get('HGNC') or
            agent.db_refs.get('CHEBI')):
        if agent.db_refs.get('TEXT'):
            name_for_node = agent.db_refs['TEXT']
            return name_for_node
    name_for_node = agent.name
    return name_for_node

def _get_node_key(agent):
    return agent.matches_key()
