import itertools
from indra.statements import *
import pygraphviz

class GraphAssembler():
    """The Graph assembler assembles INDRA Statements into a
    Graphviz node-edge graph.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be added to the assembler's list
        of Statements.

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
    """
    def __init__(self, stmts=None):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts
        self.existing_nodes = []
        self.existing_edges = []
        self.graph = pygraphviz.AGraph(directed=True, rankdir='LR',
                                       splines='spline')

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
        for stmt in self.statements:
            if isinstance(stmt, Phosphorylation):
                if stmt.enz is None:
                    continue
                self._add_node(stmt.enz)
                self._add_node(stmt.sub)
                self._add_phosphorylation(stmt.enz, stmt.sub)
            elif isinstance(stmt, Dephosphorylation):
                if stmt.enz is None:
                    continue
                self._add_node(stmt.enz)
                self._add_node(stmt.sub)
                self._add_dephosphorylation(stmt.enz, stmt.sub)
            elif isinstance(stmt, Activation):
                self._add_node(stmt.subj)
                self._add_node(stmt.obj)
                self._add_activation(stmt.subj, stmt.obj,
                                     stmt.is_activation)
            elif isinstance(stmt, Complex):
                for m in stmt.members:
                    self._add_node(m)
                self._add_complex(stmt.members)

    def get_string(self):
        """Return the assembled graph as a string.

        Returns
        -------
        graph_string : str
            The assembled graph as a string.
        """
        graph_string = self.graph.string()
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

    def save_pdf(self, file_name='graph.pdf', prog='circo'):
        """Draw the graph and save as an image or pdf file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the graph as. Default: graph.pdf
        prog : Optional[str]
            The graphviz program to use for graph layout. Default: circo
        """
        self.graph.draw(file_name, prog=prog)

    def _add_edge(self, source, target, **kwargs):
        """Add an edge to the graph."""
        style = 'solid'
        self.graph.add_edge(source, target, **kwargs)

    def _add_node(self, agent):
        """Add an Agent as a node to the graph."""
        node_name = agent.name
        if node_name in self.existing_nodes:
            return
        self.existing_nodes.append(node_name)
        node_label = agent.name
        color = "#ffeeee"
        self.graph.add_node(node_name,
                        label=node_label,
                        shape='Mrecord',
                        fillcolor=color, style="filled", color="transparent",
                        fontsize="12",
                        margin="0.06,0")

    def _add_phosphorylation(self, enz, sub):
        """Assemble a Phosphorylation statement."""
        source = enz.name
        target = sub.name
        edge_key = (source, target, 'phosphorylation')
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        params = {'color': '#000000',
                  'arrowhead': 'normal'}
        self._add_edge(source, target, **params)

    def _add_dephosphorylation(self, enz, sub):
        """Assemble a Dephosphorylation statement."""
        source = enz.name
        target = sub.name
        edge_key = (source, target, 'dephosphorylation')
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        params = {'color': '#ff0000',
                  'arrowhead': 'normal'}
        self._add_edge(source, target, **params)

    def _add_activation(self, subj, obj, rel):
        """Assemble an Activation statment."""
        source = subj.name
        target = obj.name
        edge_key = (source, target, 'activation', rel)
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        color = '#000000' if rel else '#ff0000'
        arrowhead = 'vee' if rel else 'tee'
        params = {'color': color,
                  'arrowhead': arrowhead}
        self._add_edge(source, target, **params)

    def _add_complex(self, members):
        """Assemble a Complex statement."""
        params = {'color': '#000000',
                  'arrowhead': 'none'}
        for m1, m2 in itertools.combinations(members, 2):
            edge_key = (set([m1.name, m2.name]), 'complex')
            if edge_key in self.existing_edges:
                return
            self.existing_edges.append(edge_key)
            self._add_edge(m1.name, m2.name, **params)
