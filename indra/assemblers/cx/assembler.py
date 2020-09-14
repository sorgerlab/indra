import re
import json
import time
import logging
import itertools
from ndex2.nice_cx_network import NiceCXNetwork
from collections import OrderedDict
from indra.statements import *
from indra.databases import context_client, ndex_client, get_identifiers_url, \
    url_prefixes


logger = logging.getLogger(__name__)


class NiceCxAssembler(object):
    """Assembles a Nice CX network from a set of INDRA Statements.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be assembled.
    network_name : Optional[str]
        The name of the network to be assembled. Default: indra_assembled

    Attributes
    ----------
    network : ndex2.nice_cx_network.NiceCXNetwork
        A Nice CX network object that is assembled from Statements.
    """
    def __init__(self, stmts=None, network_name=None):
        self.statements = stmts if stmts else []
        self.network = NiceCXNetwork()
        self.network.set_network_attribute('name',
                                           (network_name if network_name
                                            else 'indra_assembled'))
        self.node_keys = {}

    def make_model(self, self_loops=False, network_attributes=None):
        """Return a Nice CX network object after running assembly.

        Parameters
        ----------
        self_loops : Optional[bool]
            If False, self-loops are excluded from the network. Default: False
        network_attributes : Optional[dict]
            A dictionary containing attributes to be added to the
            assembled network.

        Returns
        -------
        ndex2.nice_cx_network.NiceCXNetwork
            The assembled Nice CX network.
        """
        for stmt in self.statements:
            agents = stmt.agent_list()
            not_none_agents = [a for a in agents if a is not None]
            if len(not_none_agents) < 2:
                continue
            for a1, a2 in itertools.combinations(not_none_agents, 2):
                a1_id = self.add_node(a1)
                a2_id = self.add_node(a2)
                if not self_loops and a1_id == a2_id:
                    continue
                edge_id = self.add_edge(a1_id, a2_id, stmt)

        prefixes = {k: v for k, v in url_prefixes.items()}
        prefixes['pubmed'] = 'https://identifiers.org/pubmed:'
        self.network.set_network_attribute('@context', json.dumps(prefixes))
        if network_attributes:
            for k, v in network_attributes.items():
                self.network.set_network_attribute(k, v, 'string')
        return self.network

    def add_node(self, agent):
        """Add an Agent to the network as a node."""
        agent_key = self.get_agent_key(agent)
        # If the node already exists
        if agent_key in self.node_keys:
            return self.node_keys[agent_key]

        # If the node doesn't exist yet
        db_ns, db_id = agent.get_grounding()
        # TODO: handle more represents name spaces
        if db_ns == 'HGNC':
            represents = 'hgnc.symbol:%s' % agent.name
        else:
            represents = None
        node_id = self.network.create_node(agent.name,
                                           node_represents=represents)
        self.node_keys[agent_key] = node_id

        # Add db_refs as aliases
        db_refs_list = ['%s:%s' % (db_name, db_id)
                        for db_name, db_id in agent.db_refs.items()
                        if db_name in url_prefixes]
        if db_refs_list:
            self.network.add_node_attribute(property_of=node_id,
                                            name='aliases',
                                            values=db_refs_list,
                                            type='list_of_string')

        # Add the type of the node, inferred from grounding
        if db_ns:
            mapped_type = db_ns_type_mappings.get(db_ns)
            if mapped_type:
                self.network.add_node_attribute(property_of=node_id,
                                                name='type',
                                                values=mapped_type,
                                                type='string')

        return node_id

    def add_edge(self, a1_id, a2_id, stmt):
        """Add a Statement to the network as an edge."""
        stmt_type = stmt.__class__.__name__
        edge_id = self.network.create_edge(a1_id, a2_id, stmt_type)
        evs = []
        for ev in stmt.evidence:
            # We skip evidences with no PMID
            if not ev.pmid:
                continue
            # We take a maximum 200 character snippet of the evidence text
            if not ev.text:
                ev_txt = 'Evidence text not available.'
            elif len(ev.text) > 200:
                ev_txt = ev.text[:200] + '...'
            else:
                ev_txt = ev.text
            # Construct a clickable PMID link with the source and evidence text
            ev_str = ('<a target="_blank" '
                      'href="https://identifiers.org/pubmed:%s">'
                      'pubmed:%s</a> (%s) %s') % (ev.pmid, ev.pmid,
                                                  ev.source_api, ev_txt)
            evs.append((ev_str, 0 if ev.text is None else 1))
        # Reorder to have ones with text first
        evs = sorted(evs, key=lambda x: x[1], reverse=True)
        # Cap at 10 pieces of evidence
        evs = [e[0] for e in evs[:10]]
        self.network.set_edge_attribute(edge_id, 'citation', evs,
                                        type='list_of_string')
        return edge_id

    def print_model(self):
        """Return the CX string of the assembled model."""
        return self.network.to_cx()

    @staticmethod
    def get_agent_key(agent):
        return agent.name


db_ns_type_mappings = {'HGNC': 'gene',
                       'UP': 'protein',
                       'FPLX': 'proteinfamily',
                       'CHEBI': 'chemical',
                       'GO': 'biological_process'}


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
            elif isinstance(stmt, Influence):
                self._add_influence(stmt)
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
                   'nodeAttributes', 'cartesianLayout']
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

    def upload_model(self, ndex_cred=None, private=True, style='default'):
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
        style : Optional[str]
            This optional parameter can either be (1)
            The UUID of an existing NDEx network whose style should be applied
            to the new network. (2) Unspecified or 'default' to use
            the default INDRA-assembled network style. (3) None to
            not set a network style.

        Returns
        -------
        network_id : str
            The UUID of the NDEx network that was created by uploading
            the assembled CX model.
        """
        cx_str = self.print_cx()
        if not ndex_cred:
            username, password = ndex_client.get_default_ndex_cred({})
            ndex_cred = {'user': username,
                         'password': password}
        network_id = ndex_client.create_network(cx_str, ndex_cred, private)
        if network_id and style:
            template_id = None if style == 'default' else style
            nretries = 3
            for retry_idx in range(nretries):
                time.sleep(3)
                try:
                    ndex_client.set_style(network_id, ndex_cred, template_id)
                    break
                except Exception:
                    msg = 'Style setting failed, '
                    if retry_idx + 1 < nretries:
                        logger.info(msg + 'retrying %d more times' %
                                    (nretries - (retry_idx+1)))
                    else:
                        logger.info(msg + 'the network will be missing style '
                                    'information.')
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
        # Here we do some bookkeeping to handle the special case where
        # a member appears twice in a complex e.g.
        # Complex(CDK12(), RECQL4(), RECQL4(), Ku())
        # and we don't want to have duplicate edges.
        added_edges = set()
        for m1, m2 in itertools.combinations(stmt.members, 2):
            m1_id = self._add_node(m1)
            m2_id = self._add_node(m2)
            if (m1_id, m2_id) not in added_edges:
                self._add_edge(m1_id, m2_id, 'Complex', stmt)
                added_edges.add((m1_id, m2_id))

    def _add_regulation(self, stmt):
        if stmt.subj is None:
            return
        subj_id = self._add_node(stmt.subj)
        obj_id = self._add_node(stmt.obj)
        stmt_type = stmt.__class__.__name__
        self._add_edge(subj_id, obj_id, stmt_type, stmt)

    def _add_influence(self, stmt):
        subj_id = self._add_node(stmt.subj.concept)
        obj_id = self._add_node(stmt.obj.concept)
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
                'i': interaction.lower()}
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
                              'n': '__INDRA json',
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
        texts = [_fix_evidence_text(e.text) for e in stmt.evidence if e.text]
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
                          'n': 'belief',
                          'v': belief_str}
        self.cx['edgeAttributes'].append(edge_attribute)

        # NOTE: supports and edgeSupports are currently
        # not shown on NDEx therefore we add text evidence as a generic
        # edgeAttribute
        if texts:
            text = texts[0]
            edge_attribute = {'po': edge_id,
                              'n': 'text',
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
    dbs = ['bel', 'biopax', 'phosphosite', 'biogrid', 'signor', 'tas', 'hprd',
           'trrust', 'ctd', 'virhostnet', 'phosphoelm', 'drugbank', 'omnipath']
    readers = ['reach', 'trips', 'sparser', 'r3', 'eidos', 'geneways', 'tees',
               'rlimsp', 'medscan']
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
    elif isinstance(stmt, Influence):
        edge_type = 'Influence'
        if stmt.overall_polarity() == -1:
            edge_polarity = 'negative'
        elif stmt.overall_polarity() == 1:
            edge_polarity = 'positive'
        else:
            edge_polarity = 'none'
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
    mir_id = agent.db_refs.get('MIRBASEM') or agent.db_refs.get('MIRBASE')
    if hgnc_id or uniprot_id:
        agent_type = 'protein'
    elif pfam_id or fa_id or be_id:
        agent_type = 'proteinfamily'
    elif chebi_id or pubchem_id:
        agent_type = 'chemical'
    elif go_id:
        agent_type = 'bioprocess'
    elif mir_id:
        agent_type = 'microrna'
    else:
        agent_type = 'other'
    return agent_type


def _fix_evidence_text(txt):
    """Eliminate some symbols to have cleaner supporting text."""
    txt = re.sub('[ ]?\( xref \)', '', txt)
    # This is to make [ xref ] become [] to match the two readers
    txt = re.sub('\[ xref \]', '[]', txt)
    txt = re.sub('[\(]?XREF_BIBR[\)]?[,]?', '', txt)
    txt = re.sub('[\(]?XREF_FIG[\)]?[,]?', '', txt)
    txt = re.sub('[\(]?XREF_SUPPLEMENT[\)]?[,]?', '', txt)
    txt = txt.strip()
    return txt


