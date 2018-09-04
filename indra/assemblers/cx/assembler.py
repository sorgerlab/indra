from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import io
import re
import json
import logging
import itertools
from collections import OrderedDict
from indra.statements import *
from indra.databases import context_client, ndex_client, get_identifiers_url

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('cx_assembler')


class CxAssembler(object):
    """This class assembles a CX network from a set of INDRA Statements.

    The CX format is an aspect oriented data mode for networks.
    The format is defined at http://www.home.ndexbio.org/data-model/.
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
    def __init__(self, stmts=None, network_name=None):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts
        if network_name is None:
            self.network_name = 'indra_assembled'
        else:
            self.network_name = network_name
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

    def make_model(self, add_indra_json=True):
        """Assemble the CX network from the collected INDRA Statements.

        This method assembles a CX network from the set of INDRA Statements.
        The assembled network is set as the assembler's cx argument.

        Parameters
        ----------
        add_indra_json : Optional[bool]
            If True, the INDRA Statement JSON annotation is added to each
            edge in the network. Default: True

        Returns
        -------
        cx_str : str
            The json serialized CX model.
        """
        self.add_indra_json = add_indra_json
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                self._add_modification(stmt)
            if isinstance(stmt, SelfModification):
                self._add_self_modification(stmt)
            elif isinstance(stmt, RegulateActivity) or \
                isinstance(stmt, RegulateAmount):
                self._add_regulation(stmt)
            elif isinstance(stmt, Complex):
                self._add_complex(stmt)
            elif isinstance(stmt, Gef):
                self._add_gef(stmt)
            elif isinstance(stmt, Gap):
                self._add_gap(stmt)
        network_description = ''
        self.cx['networkAttributes'].append({'n': 'name',
                                             'v': self.network_name})
        self.cx['networkAttributes'].append({'n': 'description',
                                             'v': network_description})
        cx_str = self.print_cx()
        return cx_str

    def print_cx(self, pretty=True):
        """Return the assembled CX network as a json string.

        Parameters
        ----------
        pretty : bool
            If True, the CX string is formatted with indentation (for human
            viewing) otherwise no indentation is used.

        Returns
        -------
        json_str : str
            A json formatted string representation of the CX network.
        """
        def _get_aspect_metadata(aspect):
            count = len(self.cx.get(aspect)) if self.cx.get(aspect) else 0
            if not count:
                return None
            data = {'name': aspect,
                    'idCounter': self._id_counter,
                    'consistencyGroup': 1,
                    'elementCount': count}
            return data
        full_cx = OrderedDict()
        full_cx['numberVerification'] = [{'longNumber': 281474976710655}]
        aspects = ['nodes', 'edges', 'supports', 'citations', 'edgeAttributes',
                   'edgeCitations', 'edgeSupports', 'networkAttributes',
                   'nodeAttributes']
        full_cx['metaData'] = []
        for aspect in aspects:
            metadata = _get_aspect_metadata(aspect)
            if metadata:
                full_cx['metaData'].append(metadata)
        for k, v in self.cx.items():
            full_cx[k] = v
        full_cx['status'] = [{'error': '', 'success': True}]
        full_cx = [{k: v} for k, v in full_cx.items()]
        if pretty:
            json_str = json.dumps(full_cx, indent=2)
        else:
            json_str = json.dumps(full_cx)
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

    def upload_model(self, ndex_cred=None, private=True):
        """Creates a new NDEx network of the assembled CX model.

        To upload the assembled CX model to NDEx, you need to have
        a registered account on NDEx (http://ndexbio.org/) and have
        the `ndex` python package installed. The uploaded network
        is private by default.

        Parameters
        ----------
        ndex_cred : Optional[dict]
            A dictionary with the following entries:
            'user': NDEx user name
            'password': NDEx password
        
        private : Optional[bool]
            Whether or not the created network will be private on NDEX.

        Returns
        -------
        network_id :  str
            The UUID of the NDEx network that was created by uploading
            the assembled CX model.
        """
        cx_str = self.print_cx()
        if not ndex_cred:
            username, password = ndex_client.get_default_ndex_cred({})
            ndex_cred = {'user': username,
                         'password': password}
        network_id = ndex_client.create_network(cx_str, ndex_cred, private)
        return network_id

    def set_context(self, cell_type):
        """Set protein expression data and mutational status as node attribute

        This method uses :py:mod:`indra.databases.context_client` to get
        protein expression levels and mutational status for a given cell type
        and set a node attribute for proteins accordingly.

        Parameters
        ----------
        cell_type : str
            Cell type name for which expression levels are queried.
            The cell type name follows the CCLE database conventions.
            Example: LOXIMVI_SKIN, BT20_BREAST
        """
        node_names = [node['n'] for node in self.cx['nodes']]
        res_expr = context_client.get_protein_expression(node_names,
                                                         [cell_type])
        res_mut = context_client.get_mutations(node_names,
                                               [cell_type])
        res_expr = res_expr.get(cell_type)
        res_mut = res_mut.get(cell_type)
        if not res_expr:
            msg = 'Could not get protein expression for %s cell type.' % \
                  cell_type
            logger.warning(msg)

        if not res_mut:
            msg = 'Could not get mutational status for %s cell type.' % \
                  cell_type
            logger.warning(msg)

        if not res_expr and not res_mut:
            return

        self.cx['networkAttributes'].append({'n': 'cellular_context',
                                             'v': cell_type})
        counter = 0
        for node in self.cx['nodes']:
            amount = res_expr.get(node['n'])
            mut = res_mut.get(node['n'])
            if amount is not None:
                node_attribute = {'po': node['@id'],
                                  'n': 'expression_amount',
                                  'v': int(amount)}
                self.cx['nodeAttributes'].append(node_attribute)
            if mut is not None:
                is_mutated = 1 if mut else 0
                node_attribute = {'po': node['@id'],
                                  'n': 'is_mutated',
                                  'v': is_mutated}
                self.cx['nodeAttributes'].append(node_attribute)
            if mut is not None or amount is not None:
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

    def _add_regulation(self, stmt):
        if stmt.subj is None:
            return
        subj_id = self._add_node(stmt.subj)
        obj_id = self._add_node(stmt.obj)
        stmt_type = stmt.__class__.__name__
        self._add_edge(subj_id, obj_id, stmt_type, stmt)

    def _add_gef(self, stmt):
        gef_id = self._add_node(stmt.gef)
        ras_id = self._add_node(stmt.ras)
        stmt_type = stmt.__class__.__name__
        self._add_edge(gef_id, ras_id, stmt_type, stmt)

    def _add_gap(self, stmt):
        gap_id = self._add_node(stmt.gap)
        ras_id = self._add_node(stmt.ras)
        stmt_type = stmt.__class__.__name__
        self._add_edge(gap_id, ras_id, stmt_type, stmt)

    def _add_node(self, agent):
        node_key = agent.name
        node_id = self._existing_nodes.get(node_key)
        if node_id is not None:
            return node_id
        node_id = self._get_new_id()
        self._existing_nodes[node_key] = node_id
        node = {'@id': node_id,
                'n': agent.name}
        self.cx['nodes'].append(node)
        self._add_node_metadata(node_id, agent)
        return node_id

    def _add_node_metadata(self, node_id, agent):
        agent_type = _get_agent_type(agent)
        node_attribute = {'po': node_id,
                          'n': 'type',
                          'v': agent_type}
        self.cx['nodeAttributes'].append(node_attribute)
        for db_name, db_ids in agent.db_refs.items():
            if not db_ids:
                logger.warning('Missing db_id for %s' % agent)
                continue
            elif isinstance(db_ids, int):
                db_id = str(db_ids)
            elif isinstance(db_ids, list):
                db_id = db_ids[0][0]
            else:
                db_id = db_ids
            url = get_identifiers_url(db_name, db_id)
            if not url:
                continue
            db_name_map = {
                'UP': 'UniProt', 'PUBCHEM': 'PubChem',
                'IP': 'InterPro', 'NXPFA': 'NextProtFamily',
                'PF': 'Pfam', 'CHEBI': 'ChEBI'}
            name = db_name_map.get(db_name)
            if not name:
                name = db_name

            node_attribute = {'po': node_id,
                              'n': name,
                              'v': url}
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

        # Add INDRA JSON
        if self.add_indra_json:
            indra_stmt_json = json.dumps(stmt.to_json())
            edge_attribute = {'po': edge_id,
                              'n': 'INDRA json',
                              'v': indra_stmt_json}
            self.cx['edgeAttributes'].append(edge_attribute)

        # Add the type of statement as the edge type
        stmt_type, stmt_polarity = _get_stmt_type(stmt)
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
            text = text.replace('XREF_BIBR', '')
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

        # Add the serialized JSON INDRA Statement
        stmt_dict = stmt.to_json()
        edge_attribute = {'po': edge_id, 'n': 'indra', 'v': stmt_dict}
        self.cx['edgeAttributes'].append(edge_attribute)

        # Add support type
        support_type = _get_support_type(stmt)
        edge_attribute = {'po': edge_id, 'n': 'supportType', 'v': support_type}
        self.cx['edgeAttributes'].append(edge_attribute)


def _get_support_type(stmt):
    dbs = ['bel', 'biopax', 'phosphosite', 'biogrid']
    readers = ['reach', 'trips', 'sparser', 'r3']
    has_db = False
    has_reading = False
    for ev in stmt.evidence:
        if ev.source_api in dbs:
            has_db = True
        if ev.source_api in readers:
            has_reading = True
    if has_db and not has_reading:
        return 'database'
    elif has_db and has_db:
        return 'database and literature'
    elif not has_db and has_reading:
        return 'literature'

def _get_stmt_type(stmt):
    if isinstance(stmt, AddModification):
        edge_type = 'Modification'
        edge_polarity = 'positive'
    elif isinstance(stmt, RemoveModification):
        edge_type = 'Modification'
        edge_polarity = 'negative'
    elif isinstance(stmt, SelfModification):
        edge_type = 'SelfModification'
        edge_polarity = 'positive'
    elif isinstance(stmt, Complex):
        edge_type = 'Complex'
        edge_polarity = 'none'
    elif isinstance(stmt, Activation):
        edge_type = 'Activation'
        edge_polarity = 'positive'
    elif isinstance(stmt, Inhibition):
        edge_type = 'Inhibition'
        edge_polarity = 'negative'
    elif isinstance(stmt, DecreaseAmount):
        edge_type = 'DecreaseAmount'
        edge_polarity = 'negative'
    elif isinstance(stmt, IncreaseAmount):
        edge_type = 'IncreaseAmount'
        edge_polarity = 'positive'
    elif isinstance(stmt, Gef):
        edge_type = 'Gef'
        edge_polarity = 'positive'
    elif isinstance(stmt, Gap):
        edge_type = 'Gap'
        edge_polarity = 'negative'
    else:
        edge_type = stmt.__class__.__str__()
        edge_polarity = 'none'
    return edge_type, edge_polarity

def _get_agent_type(agent):
    hgnc_id = agent.db_refs.get('HGNC')
    uniprot_id = agent.db_refs.get('UP')
    pfam_id = agent.db_refs.get('PF')
    fa_id = agent.db_refs.get('FA')
    chebi_id = agent.db_refs.get('CHEBI')
    pubchem_id = agent.db_refs.get('PUBCHEM')
    be_id = agent.db_refs.get('FPLX')
    go_id = agent.db_refs.get('GO')
    if hgnc_id or uniprot_id:
        agent_type = 'protein'
    elif pfam_id or fa_id or be_id:
        agent_type = 'proteinfamily'
    elif chebi_id or pubchem_id:
        agent_type = 'chemical'
    elif go_id:
        agent_type = 'bioprocess'
    else:
        agent_type = 'other'
    return agent_type
