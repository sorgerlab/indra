import networkx as nx
from indra.statements import *

class SifAssembler(object):
    def __init__(self, stmts=None):
        if stmts is None:
            self.stmts = []
        else:
            self.stmts = stmts
        self.graph = nx.DiGraph()
        self.nodes = {}

    def make_model(self):
        for st in self.stmts:
            if isinstance(st, Activation):
                s = self.add_node(st.subj)
                t = self.add_node(st.obj)
                if st.is_activation:
                    self.add_edge(s, t, {'polarity': 'positive'})
                else:
                    self.add_edge(s, t, {'polarity': 'negative'})

    def print_boolean_net(self, out_file=None):
        init_str = ''
        for node_key in self.graph.nodes():
            node_name = self.graph.node[node_key]['name']
            init_str += '%s = False\n' % node_name
        rule_str = ''
        for node_key in self.graph.nodes():
            node_name = self.graph.node[node_key]['name']
            in_edges = self.graph.in_edges(node_key)
            if not in_edges:
                continue
            parents = [e[0] for e in in_edges]
            polarities = [self.graph.edge[e[0]][node_key]['polarity']
                          for e in in_edges]
            pos_parents = [par for par, pol in zip(parents, polarities) if
                           pol == 'positive']
            neg_parents = [par for par, pol in zip(parents, polarities) if
                           pol == 'negative']

            rhs_pos_parts = []
            for par in pos_parents:
                rhs_pos_parts.append(self.graph.node[par]['name'])
            rhs_pos_str = ' or '.join(rhs_pos_parts)

            rhs_neg_parts = []
            for par in neg_parents:
                rhs_neg_parts.append(self.graph.node[par]['name'])
            rhs_neg_str = ' or '.join(rhs_neg_parts)

            if rhs_pos_str:
                if rhs_neg_str:
                    rhs_str = '(' + rhs_pos_str + \
                              ') and not (' + rhs_neg_str + ')'
                else:
                    rhs_str = rhs_pos_str
            else:
                rhs_str = 'not (' + rhs_neg_str + ')'

            node_eq = '%s* = %s\n' % (node_name, rhs_str)
            rule_str += node_eq
        full_str = init_str + '\n' + rule_str
        if out_file is not None:
            with open(out_file, 'wt') as fh:
                fh.write(full_str)
        return full_str

    def add_node(self, agent):
        node_key = agent.matches_key()
        self.graph.add_node(node_key, name=agent.name)
        return node_key

    def add_edge(self, s, t, edge_attributes=None):
        if edge_attributes is None:
            self.graph.add_edge(s, t)
        else:
            self.graph.add_edge(s, t, edge_attributes)
