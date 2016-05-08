import re
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
        self.network_name = 'indra_assembled'
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
        network_description = ''
        self.cx['networkAttributes'].append({'n': 'name',
                                             'v': self.network_name})
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
        self.add_node_metadata(node_id, agent)
        return node_id

    def add_node_metadata(self, node_id, agent):
        agent_type = get_agent_type(agent)
        node_attribute = {'po': node_id,
                          'n': 'type',
                          'v': agent_type}
        self.cx['nodeAttributes'].append(node_attribute)
        for db_name, db_ids in agent.db_refs.iteritems():
            if isinstance(db_ids, basestring):
                db_id = db_ids
            elif isinstance(db_ids, int):
                db_id = str(db_ids)
            else:
                db_id = db_ids[0]
            node_attribute = {'po': node_id,
                              'n': db_name,
                              'v': db_id}
            self.cx['nodeAttributes'].append(node_attribute)

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
        # Add the type of statement as the edge type
        stmt_type, stmt_polarity = get_stmt_type(stmt)
        edge_attribute = {'po': edge_id,
                          'n': 'type',
                          'v': stmt_type}
        self.cx['edgeAttributes'].append(edge_attribute)
        edge_attribute = {'po': edge_id,
                          'n': 'polarity',
                          'v': stmt_polarity}
        self.cx['edgeAttributes'].append(edge_attribute)
        # Add the citations for the edge
        pmids = [e.pmid for e in stmt.evidence if e.pmid]
        edge_citations = []
        for pmid in pmids:
            pmid_txt = None
            if re.match('[0-9]+', pmid):
                pmid_txt = 'pmid:' + pmid
            else:
                m = re.match('.*pubmed:([0-9]+)', pmid)
                if m:
                    pmid_txt = 'pmid:' + m.groups()[0]
                m = re.match('.*pmid:([0-9]+)', pmid)
                if m:
                    pmid_txt = 'pmid:' + m.groups()[0]
            if pmid_txt is None:
                pmid_txt = pmid
            citation_id = self._get_new_id()
            citation = {'@id': citation_id,
                        'dc:identifier': pmid_txt}
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

        # NOTE: supports and edgeSupports are currently
        # not shown on NDEx therefore we add text evidence as a generic
        # edgeAttribute
        if texts:
            text = texts[0]
            edge_attribute = {'po': edge_id,
                              'n': 'Text',
                              'v': text}
            self.cx['edgeAttributes'].append(edge_attribute)

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

def get_stmt_type(stmt):
    if isinstance(stmt, Modification):
        edge_type = 'Modification'
        edge_polarity = 'positive'
    elif isinstance(stmt, SelfModification):
        edge_type = 'SelfModification'
        edge_polarity = 'positive'
    elif isinstance(stmt, Complex):
        edge_type = 'Complex'
        edge_polarity = 'none'
    elif isinstance(stmt, ActivityActivity):
        edge_type = 'ActivityActivity'
        if stmt.relationship == 'increases':
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

def get_agent_type(agent):
    hgnc_id = agent.db_refs.get('HGNC')
    uniprot_id = agent.db_refs.get('UP')
    pfam_id = agent.db_refs.get('PF')
    fa_id = agent.db_refs.get('FA')
    chebi_id = agent.db_refs.get('CHEBI')
    if hgnc_id or uniprot_id:
        agent_type = 'protein'
    elif pfam_id or fa_id:
        agent_type = 'proteinfamily'
    elif chebi_id:
        agent_type = 'chemical'
    else:
        agent_type = 'other'
    return agent_type
