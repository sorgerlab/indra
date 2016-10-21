from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import logging
import networkx
import itertools
from indra.statements import *

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('cyjs_assembler')

class CyJSAssembler(object):
    def __init__(self, stmts=None):
        if not stmts:
            self.statements = []
        else:
            self.statements = stmts
        self._edges = []
        self._nodes = []
        self._graph = networkx.DiGraph()
        self._existing_nodes = {}
        self._id_counter = 0

    def add_statements(self, stmts):
        """Add INDRA Statements to the assembler's list of statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of :py:class:`indra.statements.Statement`
            to be added to the statement list of the assembler.
        """
        for stmt in stmts:
            self.statements.append(stmt)

    def make_model(self):
        for stmt in self.statements:
            if isinstance(stmt, Activation):
                self._add_activation(stmt)
            elif isinstance(stmt, Complex):
                self._add_complex(stmt)
            else:
                logger.warning('Unhandled statement type: %s' %
                               stmt.__class_.__name__)
        # Check for nodes that should be grouped
        #   look at self._graph to find nodes with same incoming and outgoing
        #   edges and group them
        # Create groups for those nodes
        # Merge edges for grouped nodes
        # Store in an internal data structure

    def print_cyjs(self):
        cyjs_dict = {'edges': self._edges, 'nodes': self._nodes}
        cyjs_str = json.dumps(cyjs_dict, indent=1)
        return cyjs_str

    def save_model(self, fname='model.js'):
        cyjs_str = self.print_cyjs()
        s = 'var modelElements = %s;' % cyjs_str
        with open(fname, 'wt') as fh:
            fh.write(s)

    def _add_activation(self, stmt):
        edge_type, edge_polarity = _get_stmt_type(stmt)
        edge_id = self._get_new_id()
        source_id = self._add_node(stmt.subj)
        target_id = self._add_node(stmt.obj)
        edge = {'data': {'i': edge_type, 'id': edge_id,
                         'source': source_id, 'target': target_id,
                         'polarity': edge_polarity}}
        self._edges.append(edge)
        self._graph.add_edge(source_id, target_id)

    def _add_complex(self, stmt):
        edge_type, edge_polarity = _get_stmt_type(stmt)
        for m1, m2 in itertools.combinations(stmt.members, 2):
            m1_id = self._add_node(m1)
            m2_id = self._add_node(m2)

            edge_id = self._get_new_id()
            edge = {'data': {'i': edge_type, 'id': edge_id,
                             'source': m1_id, 'target': m2_id,
                             'polarity': edge_polarity}}
            self._edges.append(edge)
            self._graph.add_edge(m1_id, m2_id)

    def _add_node(self, agent):
        node_key = agent.name
        node_id = self._existing_nodes.get(node_key)
        if node_id is not None:
            return node_id
        db_refs = _get_db_refs(agent)
        node_id = self._get_new_id()
        self._existing_nodes[node_key] = node_id
        node_name = agent.name
        node = {'data': {'id': node_id, 'name': node_name, 'db_refs': db_refs}}
        self._nodes.append(node)
        self._graph.add_node(node_id)
        return node_id

    def _get_new_id(self):
        ret = self._id_counter
        self._id_counter += 1
        return ret

    def _build_node_dict(self):
        node_ids = [x['data']['id'] for x in self._nodes]
        edge_data = [x['data'] for x in self._edges]
        node_dict = {}
        for n in node_ids:
            n_sources = []
            n_targets = []
            for e in edge_data:
                if e['target'] == n:
                    keys = ['i','polarity','source']
                    n_sources.append(_limit_dict(e, keys))
                if e['source'] == n:
                    keys = ['i','polarity','target']
                    n_targets.append(_limit_dict(e, keys))
            node_dict[n] = [n_sources,n_targets]
        return node_dict

    @staticmethod
    def _identical_nodes_lists(node_dict):
        identicals = []
        for n1 in node_dict:
            for n2 in node_dict:
                if n1 != n2:
                    if node_dict[n1][0] == node_dict[n2][0]:
                        if node_dict[n1][1] == node_dict[n2][1]:
                            identicals.append([n1,n2])
        return identicals

    @staticmethod
    def _combine_identicals(identicals):

        combined = list(sorted(identicals))
        sets = []
        for a in combined:
            combos = []
            trigger1 = False
            for b in combined:
                if a!=b:
                    c = list(set(a+b))
                    if len(c) < (len(a+b)):
                        trigger1 = True
                        combos = list(sorted(list(set(combos+c))))
                    if len([x for x in a if x in b]) > 0:
                        combos = list(sorted(list(set(combos+c))))
                        trigger1 = True
            if trigger1 == False:
                combos = list(sorted(list(set(combos + a))))
            sets.append(combos)
        node_groupings = []
        for x in sets:
            if x not in node_groupings:
                node_groupings.append(x)
        return node_groupings

    def node_sets(self):
        node_dict = self._build_node_dict()
        identicals = self._identical_nodes_lists(node_dict)
        node_groupings = list(sorted(identicals))

        while True:
            pre_len = len(node_groupings)
            node_groupings = self._combine_identicals(node_groupings)
            post_len = len(node_groupings)

            if pre_len == post_len:
                break

        return node_groupings

    def _collect_nodes(self):


        new_edges = []
        group_nodes = []

        for g in node_groupings:

            # make new group node

            new_group_node = {'data': {'id': (self._get_new_id()),
                                       'name': ('Group'+str(g)),
                                       'parent': ''}}

            self._nodes.append(new_group_node)

            # kill old edges, get their info
            # make new edges to group node with old info

            kill_edges = []
            kill_nodes = []

            source_list = []
            target_list = []

            for n in g:
                edges = self._edges

                for i in range(0, len(self._nodes)):
                    if self._nodes[i]['data']['id'] == n:
                        self._nodes[i]['data']['parent'] = 'Group'+str(g)


                for e in edges:
                    if e['data']['target'] == n:
                        if e not in kill_edges:
                            kill_edges.append(e)
                        if e['data']['source'] not in source_list:
                            source_list.append(e['data']['source'])
                            new_source_edge = e
                            new_source_edge['data']['target'] = new_group_node['data']['id']
                            new_source_edge['data']['id'] = self._get_new_id()
                            new_edges.append(new_source_edge)



                    if e['data']['source'] == n:
                        if e not in kill_edges:
                            kill_edges.append(e)
                        if e['data']['target'] not in target_list:
                            target_list.append(e['data']['target'])
                            new_target_edge = e
                            new_target_edge['data']['source'] = new_group_node['data']['id']
                            new_target_edge['data']['id'] = self._get_new_id()
                            new_edges.append(new_target_edge)

            for e in kill_edges:
                self._edges.remove(e)

            for e in new_edges:
                self._edges.append(e)

        return None


def _get_db_refs(agent):
    cyjs_db_refs = {}
    for db_name, db_ids in agent.db_refs.items():
        if isinstance(db_ids, int):
            db_id = str(db_ids)
        elif isinstance(db_ids, basestring):
            db_id = db_ids
        else:
            db_id = db_ids[0]
        if db_name == 'UP':
            name = 'UniProt'
            val = 'http://identifiers.org/uniprot/%s' % db_id
        elif db_name == 'HGNC':
            name = 'HGNC'
            val = 'http://identifiers.org/hgnc/HGNC:%s' % db_id
        elif db_name == 'CHEBI':
            name = 'ChEBI'
            val = 'http://identifiers.org/chebi/%s' % db_id
        elif db_name == 'PUBCHEM':
            name = 'PubChem'
            val = 'http://identifiers.org/pubchem.compound/%s' % db_id
        elif db_name == 'HMDB':
            name = 'HMDB'
            val = 'http://identifiers.org/hmdb/%s' % db_id
        elif db_name == 'GO':
            name = 'GO'
            val = 'http://identifiers.org/go/%s' % db_id
        elif db_name == 'MESH':
            name = 'MESH'
            val = 'http://identifiers.org/mesh/%s' % db_id
        elif db_name == 'IP':
            name = 'InterPro'
            val = 'http://identifiers.org/interpro/%s' % db_id
        elif db_name == 'TEXT':
            continue
        else:
            val = db_id
            name = db_name
        cyjs_db_refs[name] = val
    return cyjs_db_refs

def _get_stmt_type(stmt):
    if isinstance(stmt, Modification):
        edge_type = 'Modification'
        edge_polarity = 'positive'
    elif isinstance(stmt, SelfModification):
        edge_type = 'SelfModification'
        edge_polarity = 'positive'
    elif isinstance(stmt, Complex):
        edge_type = 'Complex'
        edge_polarity = 'none'
    elif isinstance(stmt, Activation):
        edge_type = 'Activation'
        if stmt.is_activation:
            edge_polarity = 'positive'
        else:
            edge_polarity = 'negative'
    elif isinstance(stmt, RasGef):
        edge_type = 'RasGef'
        edge_polarity = 'positive'
    elif isinstance(stmt, RasGap):
        edge_type = 'RasGap'
        edge_polarity = 'negative'
    else:
        edge_type = stmt.__class__.__str__()
        edge_polarity = 'none'
    return edge_type, edge_polarity

def _limit_dict(d, keys):
    d2 = { k: d[k] for k in keys }
    return d2
