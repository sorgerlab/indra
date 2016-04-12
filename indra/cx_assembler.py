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
                   'citations': [], 'edgeCitations': [],
                   'supports': [], 'edgeSupports': [],
                   'networkAttributes': []}
        self.id_counter = 0

    def _get_new_id(self):
        ret = self.id_counter
        self.id_counter += 1
        return ret

    def add_statements(self, stmts):
        for stmt in stmts:
            self.statements.append(stmt)

    def make_model(self):
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                self.add_modification(stmt)
            if isinstance(stmt, SelfModification):
                self.add_self_modification(stmt)
            elif isinstance(stmt, ActivityActivity):
                self.add_activityactivity(stmt)
            elif isinstance(stmt, Complex):
                self.add_complex(stmt)
            elif isinstance(stmt, RasGef):
                self.add_rasgef(stmt)
            elif isinstance(stmt, RasGap):
                self.add_rasgap(stmt)
        network_name = 'indra_assembled'
        network_description = ''
        self.cx['networkAttributes'].append({'n': 'name', 'v': network_name})
        self.cx['networkAttributes'].append({'n': 'description',
                                             'v': network_description})

    def add_modification(self, stmt):
        if stmt.enz is None:
            return
        enz_id = self.add_node(stmt.enz)
        sub_id = self.add_node(stmt.sub)
        stmt_type = stmt.__class__.__name__
        self.add_edge(enz_id, sub_id, stmt_type, stmt)

    def add_self_modification(self, stmt):
        enz_id = self.add_node(stmt.enz)
        stmt_type = stmt.__class__.__name__
        self.add_edge(enz_id, enz_id, stmt_type, stmt)

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

    def add_rasgef(self, stmt):
        gef_id = self.add_node(stmt.gef)
        ras_id = self.add_node(stmt.ras)
        stmt_type = stmt.__class__.__name__
        self.add_edge(gef_id, ras_id, stmt_type, stmt)

    def add_rasgap(self, stmt):
        gap_id = self.add_node(stmt.gap)
        ras_id = self.add_node(stmt.ras)
        stmt_type = stmt.__class__.__name__
        self.add_edge(gap_id, ras_id, stmt_type, stmt)

    def add_node(self, agent):
        node_key = agent.name
        try:
            node_id = self.existing_nodes[node_key]
            return node_id
        except KeyError:
            pass
        node_id = self._get_new_id()
        self.existing_nodes[node_key] = node_id
        node = {'@id': node_id,
                'n': agent.name}
        self.cx['nodes'].append(node)
        for db_name, db_ids in agent.db_refs.iteritems():
            node_attribute = {'po': node_id,
                              'n': db_name,
                              'v': db_ids}
            self.cx['nodeAttributes'].append(node_attribute)
        return node_id

    def add_edge(self, source, target, interaction, stmt):
        edge_key = (source, target, interaction)
        try:
            edge_id = self.existing_edges[edge_key]
            return edge_id
        except KeyError:
            pass
        edge_id = self._get_new_id()
        self.existing_nodes[edge_key] = edge_id
        edge = {'@id': edge_id,
                's': source,
                't': target,
                'i': interaction}
        self.cx['edges'].append(edge)
        self.add_edge_metadata(edge_id, stmt)
        return edge_id

    def add_edge_metadata(self, edge_id, stmt):
        '''
        Add all annotations, evidence, citations, etc. to the edge
        '''
        # Add the string of the statement itself
        indra_stmt_str = '%s' % stmt
        edge_attribute = {'po': edge_id,
                          'n': 'INDRA statement',
                          'v': indra_stmt_str}
        self.cx['edgeAttributes'].append(edge_attribute)

        # Add the citations for the edge
        pmids = [e.pmid for e in stmt.evidence if e.pmid]
        edge_citations = []
        for pmid in pmids:
            citation_id = self._get_new_id()
            citation = {'@id': citation_id,
                        'dc:identifier': 'pmid:%s' % pmid}
            self.cx['citations'].append(citation)
            edge_citations.append(citation_id)
        if edge_citations:
            edge_citation = {'citations': edge_citations,
                             'po': [edge_id]}
            self.cx['edgeCitations'].append(edge_citation)

        # Add the textual supports for the edge
        texts = [e.text for e in stmt.evidence if e.text]
        edge_supports = []
        for text in texts:
            support_id = self._get_new_id()
            support = {'@id': support_id,
                       'text': text}
            self.cx['supports'].append(support)
            edge_supports.append(support_id)
        if edge_supports:
            edge_support = {'supports': edge_supports,
                            'po': [edge_id]}
            self.cx['edgeSupports'].append(edge_support)


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
