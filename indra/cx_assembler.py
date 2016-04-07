import json
import itertools
from indra.statements import *

class CxAssembler():
    # http://www.ndexbio.org/data-model/
    def __init__(self):
        self.statements = []
        self.existing_nodes = []
        self.existing_edges = []
        self.cx = {'nodes': [], 'edges': [],
                      'nodeAttributes': [], 'edgeAttributes': []}
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

    def add_phosphorylation(self, stmt):
        if stmt.enz is None:
            return
        enz_id = self.add_node(stmt.enz)
        sub_id = self.add_node(stmt.sub)
        self.add_edge(enz_id, sub_id, 'Phosphorylation', stmt)

    def add_node(self, agent):
        node_id = self.id_counter
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
        edge_id = self.id_counter
        edge = {'@id': edge_id,
                's': source,
                't': target,
                'i': interaction}
        self.cx['edges'].append(edge)
        if stmt.evidence:
            if stmt.evidence[0].pmid:
                edge_attribute = {'po': edge_id,
                                  'n': 'pmid',
                                  'v': stmt.evidence[0].pmid}
                self.cx['edgeAttributes'].append(edge_attribute)
        return edge_id

    def print_json(self):
        json_str = json.dumps(self.cx)
        return json_str
