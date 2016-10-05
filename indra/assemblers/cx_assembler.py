from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import json
import itertools
from collections import OrderedDict
from indra.statements import *
from indra.databases import context_client

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

class CxAssembler():
    """This class assembles a CX network from a set of INDRA Statements.

    The CX format is an aspect oriented data mode for networks.
    The format is defined at http://www.ndexbio.org/data-model/.
    The CX format is the standard for NDEx and is compatible with
    CytoScape via the CyNDEx plugin.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be assembled.
    network_name : Optional[str]
        The name of the network to be assembled. Default: indra_assembled

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to be assembled.
    network_name : str
        The name of the network to be assembled.
    cx : dict
        The structure of the CX network that is assembled.
    """
    def __init__(self, stmts=None, network_name='indra_assembled'):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts
        self.network_name = 'indra_assembled'
        self.cx = {'nodes': [], 'edges': [],
                   'nodeAttributes': [], 'edgeAttributes': [],
                   'citations': [], 'edgeCitations': [],
                   'supports': [], 'edgeSupports': [],
                   'networkAttributes': []}
        self._existing_nodes = {}
        self._existing_edges = {}
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
        """Assemble the CX network from the collected INDRA Statements.

        This method assembles a CX network from the set of INDRA Statements.
        The assembled network is set as the assembler's cx argument.
        """
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                self._add_modification(stmt)
            if isinstance(stmt, SelfModification):
                self._add_self_modification(stmt)
            elif isinstance(stmt, Activation):
                self._add_activation(stmt)
            elif isinstance(stmt, Complex):
                self._add_complex(stmt)
            elif isinstance(stmt, RasGef):
                self._add_rasgef(stmt)
            elif isinstance(stmt, RasGap):
                self._add_rasgap(stmt)
        network_description = ''
        self.cx['networkAttributes'].append({'n': 'name',
                                             'v': self.network_name})
        self.cx['networkAttributes'].append({'n': 'description',
                                             'v': network_description})
    def print_cx(self):
        """Return the assembled CX network as a json string."""
        full_cx = OrderedDict()
        full_cx['numberVerification'] = [{'longNumber': 281474976710655}]
        full_cx['metaData'] = [{'idCounter': self._id_counter,
                                'name': 'nodes'},
                               {'idCounter': self._id_counter,
                                'name': 'edges'}]
        for k, v in self.cx.items():
            full_cx[k] = v
        full_cx = [{k: v} for k, v in full_cx.items()]
        json_str = json.dumps(full_cx, indent=2)
        return json_str

    def save_model(self, file_name='model.cx'):
        """Save the assembled CX network in a file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the CX network to. Default: model.cx
        """
        with open(file_name, 'wt') as fh:
            cx_str = self.print_cx()
            fh.write(cx_str)

    def set_context(self, cell_type):
        """Set protein expression data as node attribute

        This method uses :py:mod:`indra.databases.context_client` to get
        protein expression levels for a given cell type and set a node
        attribute for proteins accordingly.

        Parameters
        ----------
        cell_type : str
            Cell type name for which expression levels are queried.
            The cell type name follows the CCLE database conventions.

        Example: LOXIMVI_SKIN, BT20_BREAST
        """
        node_names = [node['n'] for node in self.cx['nodes']]
        res = context_client.get_protein_expression(node_names, cell_type)
        if not res:
            logger.warning('Could not get context for %s cell type.' %
                           cell_type)
            return
        self.cx['networkAttributes'].append({'n': 'cellular_context',
                                             'v': cell_type})
        counter = 0
        for node in self.cx['nodes']:
            amount = res.get(node['n'])
            if amount is None:
                continue
            node_attribute = {'po': node['@id'],
                              'n': 'expression_amount',
                              'v': int(amount[cell_type])}
            self.cx['nodeAttributes'].append(node_attribute)
            counter += 1
        logger.info('Set context for %d nodes.' % counter)

    def _get_new_id(self):
        ret = self._id_counter
        self._id_counter += 1
        return ret

    def _add_modification(self, stmt):
        if stmt.enz is None:
            return
        enz_id = self._add_node(stmt.enz)
        sub_id = self._add_node(stmt.sub)
        stmt_type = stmt.__class__.__name__
        self._add_edge(enz_id, sub_id, stmt_type, stmt)

    def _add_self_modification(self, stmt):
        enz_id = self._add_node(stmt.enz)
        stmt_type = stmt.__class__.__name__
        self._add_edge(enz_id, enz_id, stmt_type, stmt)

    def _add_complex(self, stmt):
        for m1, m2 in itertools.combinations(stmt.members, 2):
            m1_id = self._add_node(m1)
            m2_id = self._add_node(m2)
            self._add_edge(m1_id, m2_id, 'Complex', stmt)

    def _add_activation(self, stmt):
        subj_id = self._add_node(stmt.subj)
        obj_id = self._add_node(stmt.obj)
        # TODO: take into account relation here
        self._add_edge(subj_id, obj_id, 'Activation', stmt)

    def _add_rasgef(self, stmt):
        gef_id = self._add_node(stmt.gef)
        ras_id = self._add_node(stmt.ras)
        stmt_type = stmt.__class__.__name__
        self._add_edge(gef_id, ras_id, stmt_type, stmt)

    def _add_rasgap(self, stmt):
        gap_id = self._add_node(stmt.gap)
        ras_id = self._add_node(stmt.ras)
        stmt_type = stmt.__class__.__name__
        self._add_edge(gap_id, ras_id, stmt_type, stmt)

    def _add_node(self, agent):
        node_key = agent.name
        try:
            node_id = self._existing_nodes[node_key]
            return node_id
        except KeyError:
            pass
        node_id = self._get_new_id()
        self._existing_nodes[node_key] = node_id
        node = {'@id': node_id,
                'n': agent.name}
        self.cx['nodes'].append(node)
        self._add_node_metadata(node_id, agent)
        return node_id

    def _add_node_metadata(self, node_id, agent):
        agent_type = get_agent_type(agent)
        node_attribute = {'po': node_id,
                          'n': 'type',
                          'v': agent_type}
        self.cx['nodeAttributes'].append(node_attribute)
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

            node_attribute = {'po': node_id,
                              'n': name,
                              'v': val}
            self.cx['nodeAttributes'].append(node_attribute)

    def _add_edge(self, source, target, interaction, stmt):
        edge_key = (source, target, interaction)
        try:
            edge_id = self._existing_edges[edge_key]
            return edge_id
        except KeyError:
            pass
        edge_id = self._get_new_id()
        self._existing_nodes[edge_key] = edge_id
        edge = {'@id': edge_id,
                's': source,
                't': target,
                'i': interaction}
        self.cx['edges'].append(edge)
        self._add_edge_metadata(edge_id, stmt)
        return edge_id

    def _add_edge_metadata(self, edge_id, stmt):
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
        pmids_added = []
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
            if pmid_txt not in pmids_added:
                citation_id = self._get_new_id()
                citation = {'@id': citation_id,
                            'dc:identifier': pmid_txt}
                self.cx['citations'].append(citation)
                edge_citations.append(citation_id)
                pmids_added.append(pmid_txt)
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

        belief_str = '%.2f' % stmt.belief
        edge_attribute = {'po': edge_id,
                          'n': 'Belief score',
                          'v': belief_str}
        self.cx['edgeAttributes'].append(edge_attribute)

        # NOTE: supports and edgeSupports are currently
        # not shown on NDEx therefore we add text evidence as a generic
        # edgeAttribute
        if texts:
            text = texts[0]
            edge_attribute = {'po': edge_id,
                              'n': 'Text',
                              'v': text}
            self.cx['edgeAttributes'].append(edge_attribute)

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

def get_agent_type(agent):
    hgnc_id = agent.db_refs.get('HGNC')
    uniprot_id = agent.db_refs.get('UP')
    pfam_id = agent.db_refs.get('PF')
    fa_id = agent.db_refs.get('FA')
    chebi_id = agent.db_refs.get('CHEBI')
    be_id = agent.db_refs.get('BE')
    go_id = agent.db_refs.get('GO')
    if hgnc_id or uniprot_id:
        agent_type = 'protein'
    elif pfam_id or fa_id:
        agent_type = 'proteinfamily'
    elif chebi_id:
        agent_type = 'chemical'
    elif be_id:
        agent_type = 'proteinfamily'
    elif go_id:
        agent_type = 'bioprocess'
    else:
        agent_type = 'other'
    return agent_type
