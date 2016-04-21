import itertools
from indra.statements import *
import pygraphviz

class GraphAssembler():
    def __init__(self):
        self.statements = []
        self.existing_nodes = []
        self.existing_edges = []
        self.graph = pygraphviz.AGraph(directed=True, rankdir='LR',
                                       splines='spline')

    def add_statements(self, stmts):
        for stmt in stmts:
            self.statements.append(stmt)

    def make_model(self):
        for stmt in self.statements:
            if isinstance(stmt, Phosphorylation):
                if stmt.enz is None:
                    continue
                self.add_node(stmt.enz)
                self.add_node(stmt.sub)
                self.add_phosphorylation(stmt.enz, stmt.sub)
            elif isinstance(stmt, Dephosphorylation):
                if stmt.enz is None:
                    continue
                self.add_node(stmt.enz)
                self.add_node(stmt.sub)
                self.add_dephosphorylation(stmt.enz, stmt.sub)
            elif isinstance(stmt, ActivityActivity):
                self.add_node(stmt.subj)
                self.add_node(stmt.obj)
                self.add_activityactivity(stmt.subj, stmt.obj,
                                          stmt.relationship)
            elif isinstance(stmt, Complex):
                for m in stmt.members:
                    self.add_node(m)
                self.add_complex(stmt.members)

    def add_edge(self, source, target, **kwargs):
        style = 'solid'
        self.graph.add_edge(source, target, **kwargs)

    def add_phosphorylation(self, enz, sub):
        source = enz.name
        target = sub.name
        edge_key = (source, target, 'phosphorylation')
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        params = {'color': '#000000',
                  'arrowhead': 'normal'}
        self.add_edge(source, target, **params)

    def add_dephosphorylation(self, enz, sub):
        source = enz.name
        target = sub.name
        edge_key = (source, target, 'dephosphorylation')
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        params = {'color': '#ff0000',
                  'arrowhead': 'normal'}
        self.add_edge(source, target, **params)

    def add_activityactivity(self, subj, obj, rel):
        source = subj.name
        target = obj.name
        edge_key = (source, target, 'activity', rel)
        if edge_key in self.existing_edges:
            return
        self.existing_edges.append(edge_key)
        color = '#000000' if rel == 'increases' else '#ff0000'
        arrowhead = 'vee' if rel == 'increases' else 'tee'
        params = {'color': color,
                  'arrowhead': arrowhead}
        self.add_edge(source, target, **params)

    def add_complex(self, members):
        params = {'color': '#000000',
                  'arrowhead': 'none'}
        for m1, m2 in itertools.combinations(members, 2):
            edge_key = (set([m1.name, m2.name]), 'complex')
            if edge_key in self.existing_edges:
                return
            self.existing_edges.append(edge_key)
            self.add_edge(m1.name, m2.name, **params)

    def add_node(self, agent):
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

    def get_string(self):
        graph_string = self.graph.string()
        graph_string = graph_string.replace('\\N', '\\n')
        return graph_string

    def save_dot(self, fname='graph.dot'):
        s = self.get_string()
        with open(fname, 'wt') as fh:
            fh.write(s)

    def save_pdf(self, fname='graph.pdf'):
        self.graph.draw(fname, prog='circo')
