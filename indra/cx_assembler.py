import json
import itertools
from collections import OrderedDict
from indra.statements import *

class CxAssembler():
    # http://www.ndexbio.org/data-model/
    def __init__(self):
        self.statements = []
        self.existing_nodes = {}
        self.existing_edges = {}
        self.cx = {'nodes': [], 'edges': [],
                   'nodeAttributes': [], 'edgeAttributes': [],
                   'networkAttributes': []}
        self.id_counter = 0

    def add_statements(self, stmts):
        for stmt in stmts:
            self.statements.append(stmt)

    def make_model(self):
        for stmt in self.statements:
            if isinstance(stmt, Phosphorylation):
                self.add_phosphorylation(stmt)
            elif isinstance(stmt, Dephosphorylation):
                self.add_dephosphorylation(stmt)
            elif isinstance(stmt, ActivityActivity):
                self.add_activityactivity(stmt)
            elif isinstance(stmt, Complex):
                self.add_complex(stmt)
        network_name = 'indra_assembled'
        network_description = ''
        self.cx['networkAttributes'].append({'n': 'name', 'v': network_name})
        self.cx['networkAttributes'].append({'n': 'description',
                                             'v': network_description})

    def add_phosphorylation(self, stmt):
        if stmt.enz is None:
            return
        enz_id = self.add_node(stmt.enz)
        sub_id = self.add_node(stmt.sub)
        self.add_edge(enz_id, sub_id, 'Phosphorylation', stmt)

    def add_dephosphorylation(self, stmt):
        if stmt.enz is None:
            return
        enz_id = self.add_node(stmt.enz)
        sub_id = self.add_node(stmt.sub)
        self.add_edge(enz_id, sub_id, 'Dephosphorylation', stmt)

    def add_complex(self, stmt):
        for m1, m2 in itertools.combinations(stmt.members, 2):
            m1_id = self.add_node(m1)
            m2_id = self.add_node(m2)
            self.add_edge(m1_id, m2_id, 'Complex', stmt)

    def add_activityactivity(self, stmt):
        subj_id = self.add_node(stmt.subj)
        obj_id = self.add_node(stmt.obj)
        # TODO: take into account relation here
        self.add_edge(subj_id, obj_id, 'ActivityActivity', stmt)

    def add_node(self, agent):
        node_key = agent.name
        try:
            node_id = self.existing_nodes[node_key]
            return node_id
        except KeyError:
            pass
        node_id = self.id_counter
        self.existing_nodes[node_key] = node_id
        node = {'@id': node_id,
                'n': agent.name}
        self.cx['nodes'].append(node)
        for db_name, db_ids in agent.db_refs.iteritems():
            node_attribute = {'po': node_id,
                              'n': db_name,
                              'v': db_ids}
            self.cx['nodeAttributes'].append(node_attribute)
        self.id_counter += 1
        return node_id

    def add_edge(self, source, target, interaction, stmt):
        edge_key = (source, target, interaction)
        try:
            edge_id = self.existing_edges[edge_key]
            return edge_id
        except KeyError:
            pass
        edge_id = self.id_counter
        self.existing_nodes[edge_key] = edge_id
        edge = {'@id': edge_id,
                's': source,
                't': target,
                'i': interaction}
        self.cx['edges'].append(edge)
        indra_stmt_str = '%s' % stmt
        edge_attribute = {'po': edge_id,
                          'n': 'INDRA statement',
                          'v': indra_stmt_str}
        self.cx['edgeAttributes'].append(edge_attribute)
        self.id_counter += 1
        return edge_id

    def print_cx(self):
        full_cx = OrderedDict()
        full_cx['numberVerification'] = [{'longNumber': 281474976710655}]
        full_cx['metaData'] = [{'idCounter': self.id_counter, 'name': 'nodes'},
                               {'idCounter': self.id_counter, 'name': 'edges'}]
        for k, v in self.cx.iteritems():
            full_cx[k] = v
        full_cx = [{k: v} for k, v in full_cx.iteritems()]
        json_str = json.dumps(full_cx, indent=2)
        return json_str

    def save_model(self, fname='model.cx'):
        with open(fname, 'wt') as fh:
            cx_str = self.print_cx()
            fh.write(cx_str)
