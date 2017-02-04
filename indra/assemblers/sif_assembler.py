from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import numpy
import logging
import itertools
import networkx as nx
from indra.statements import *

logger = logging.getLogger('sif_assembler')

class SifAssembler(object):
    """The SIF assembler assembles INDRA Statements into a networkx graph.

    This graph can then be exported into SIF (simple ineraction format) or
    a Boolean network.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be added to the assembler's list
        of Statements.

    Attributes
    ----------
    graph : networkx.DiGraph
        A networkx graph that is assembled by this assembler.
    """
    def __init__(self, stmts=None):
        if stmts is None:
            self.stmts = []
        else:
            self.stmts = stmts
        self.graph = nx.DiGraph()
        self._use_name_as_key = False

    def make_model(self, use_name_as_key=False, include_mods=False,
                   include_complexes=False):
        """Assemble the graph from the assembler's list of INDRA Statements.

        Parameters
        ----------
        use_name_as_key : boolean
            If True, uses the name of the agent as the key to the nodes in
            the network. If False (default) uses the matches_key() of the
            agent.
        include_mods : boolean
            If True, adds Modification statements into the graph as directed
            edges. Default is False.
        include_complexes : boolean
            If True, creates two edges (in both directions) between all pairs
            of nodes in Complex statements. Default is False.
        """
        self.graph = nx.DiGraph()
        self._use_name_as_key = use_name_as_key
        for st in self.stmts:
            support_all = len(st.evidence)
            support_pmid = len(set([ev.pmid for ev in st.evidence
                                    if ev.pmid is not None]))
            attr = {'polarity': 'unknown', 'support_all': support_all,
                    'support_pmid': support_pmid}
            if isinstance(st, RegulateActivity):
                attr['polarity'] = ('positive' if st.is_activation
                                    else 'negative')
                self._add_node_edge(st.subj, st.obj, attr)
            elif include_mods and isinstance(st, Modification):
                self._add_node_edge(st.agent_list()[0], st.agent_list()[1], attr)
            elif include_mods and \
                 (isinstance(st, RasGap) or isinstance(st, DecreaseAmount)):
                attr['polarity'] = 'negative'
                self._add_node_edge(st.agent_list()[0], st.agent_list()[1], attr)
            elif include_mods and \
                 (isinstance(st, RasGef) or isinstance(st, IncreaseAmount)):
                attr['polarity'] = 'positive'
                self._add_node_edge(st.agent_list()[0], st.agent_list()[1], attr)
            elif include_complexes and isinstance(st, Complex):
                # Create s->t edges between all possible pairs of complex
                # members
                for node1, node2 in itertools.permutations(st.members, 2):
                    self._add_node_edge(node1, node2, attr)

    def _add_node_edge(self, s, t, attributes):
        if s is not None:
            s = self._add_node(s)
            t = self._add_node(t)
            self._add_edge(s, t, attributes)

    def set_edge_weights(self, attribute):
        max_val = 0
        for s, t, attr in self.graph.edges(data=True):
            max_val = max(attr[attribute], max_val)
        if max_val > 0:
            for s, t, attr in self.graph.edges(data=True):
                self.graph[s][t]['weight'] = \
                    1.0 / attr[attribute]

    def print_model(self):
        """Return a SIF string of the assembled model."""
        sif_str = ''
        for edge in self.graph.edges_iter(data=True):
            n1 = edge[0]
            n2 = edge[1]
            data = edge[2]
            polarity = data.get('polarity')
            if polarity == 'negative':
                rel = '-1'
            elif polarity == 'positive':
                rel = '1'
            else:
                continue
            sif_str += '%s %s %s\n' % (n1, rel, n2)
        return sif_str

    def print_boolean_net(self, out_file=None):
        """Return a Boolean network from the assembled graph.

        See https://github.com/ialbert/booleannet for details about
        the format used to encode the Boolean rules.

        Parameters
        ----------
        out_file : Optional[str]
            A file name in which the Boolean network is saved.

        Returns
        -------
        full_str : str
            The string representing the Boolean network.
        """
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

    def _add_node(self, agent):
        if self._use_name_as_key:
            node_key = agent.name
        else:
            node_key = agent.matches_key()
        self.graph.add_node(node_key, name=agent.name)
        return node_key

    def _add_edge(self, s, t, edge_attributes=None):
        if edge_attributes is None:
            self.graph.add_edge(s, t)
        else:
            self.graph.add_edge(s, t, edge_attributes)
