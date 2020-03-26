from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import logging
from indra.databases import uniprot_client, hgnc_client
from indra.statements import *


logger = logging.getLogger(__name__)


_stmt_map = {
    'phosphorylates': Phosphorylation,
    'controls-phosphorylation-of': Phosphorylation,
    'in-complex-with': Complex,
    'NEGATIVE_INFLUENCE': Inhibition,
    'POSITIVE_INFLUENCE': Activation,
    'binds': Complex,
    'interacts with': Complex,
}


def _get_dict_from_list(dict_key, list_of_dicts):
    """Retrieve a specific dict from a list of dicts.

    Parameters
    ----------
    dict_key : str
        The (single) key of the dict to be retrieved from the list.
    list_of_dicts : list
        The list of dicts to search for the specific dict.

    Returns
    -------
    dict value
        The value associated with the dict_key (e.g., a list of nodes or
        edges).
    """
    the_dict = [cur_dict for cur_dict in list_of_dicts
                if cur_dict.get(dict_key)]
    if not the_dict:
        raise ValueError('Could not find a dict with key %s' % dict_key)
    return the_dict[0][dict_key]


class NdexCxProcessor(object):
    """The NdexCxProcessor extracts INDRA Statements from Cytoscape CX JSON.

    Parameters
    ----------
    cx : list of dicts
        JSON content containing the Cytoscape network in CX format.
    summary : Optional[dict]
        The network summary object which can be obtained via
        get_network_summary through the web service. THis contains metadata
        such as the owner and the creation time of the network.

    Attributes
    ----------
    statements : list
        A list of extracted INDRA Statements. Not all edges in the network
        may be converted into Statements.
    """
    def __init__(self, cx, summary=None, require_grounding=True):
        self.cx = cx
        self.statements = []
        self.require_grounding = require_grounding
        # Initialize the dict mapping node IDs to gene names
        self._node_names = {}
        self._node_agents = {}
        self._network_info = {}
        self._edge_attributes = {}
        summary = summary if summary else {}
        self._initialize_node_attributes()
        self._initialize_node_agents()
        self._initialize_network_info(summary)
        self._initialize_edge_attributes()

    def _initialize_node_agents(self):
        """Initialize internal dicts containing node information."""
        nodes = _get_dict_from_list('nodes', self.cx)
        invalid_genes = []
        for node in nodes:
            id = node['@id']
            cx_db_refs = self.get_aliases(node)
            node_name = node['n']
            up_id = cx_db_refs.get('UP')
            if up_id:
                db_refs = {'UP': up_id, 'TEXT': node_name}
                hgnc_id = uniprot_client.get_hgnc_id(up_id)
                if hgnc_id:
                    db_refs['HGNC'] = hgnc_id
                    gene_name = hgnc_client.get_hgnc_name(hgnc_id)
                else:
                    gene_name = uniprot_client.get_gene_name(up_id)
                agent = Agent(gene_name, db_refs=db_refs)
                self._node_names[id] = gene_name
                self._node_agents[id] = agent
                continue
            else:
                self._node_names[id] = node_name
                hgnc_id = hgnc_client.get_hgnc_id(node_name)
                db_refs = {'TEXT': node_name}
                if not hgnc_id:
                    if not self.require_grounding:
                        self._node_agents[id] = \
                                Agent(node_name, db_refs=db_refs)
                    invalid_genes.append(node_name)
                else:
                    db_refs.update({'HGNC': hgnc_id})
                    up_id = hgnc_client.get_uniprot_id(hgnc_id)
                    # It's possible that a valid HGNC ID will not have a
                    # Uniprot ID, as in the case of HOTAIR (HOX transcript
                    # antisense RNA, HGNC:33510)
                    if up_id:
                        db_refs.update({'UP': up_id})
                    self._node_agents[id] = Agent(node_name, db_refs=db_refs)
        if invalid_genes:
            verb = 'Skipped' if self.require_grounding else 'Included'
            logger.info('%s invalid gene symbols: %s' %
                        (verb, ', '.join(invalid_genes)))

    def _initialize_network_info(self, summary):
        self._network_info['externalId'] = summary.get('externalId')
        self._network_info['owner'] = summary.get('owner')

    def _initialize_edge_attributes(self):
        edge_attr = _get_dict_from_list('edgeAttributes', self.cx)
        for ea in edge_attr:
            edge_id = ea.get('po')
            ea_type = ea.get('n')
            ea_value = ea.get('v')
            ea_info = self._edge_attributes.get(edge_id)
            # If we don't have any info about this edge, initialize an empty
            # dict
            if ea_info is None:
                ea_info = {'pmids': []}
                self._edge_attributes[edge_id] = ea_info
            # Collect PMIDs from the various edge types
            if ea_type == 'ndex:citation' or ea_type == 'citation_ids':
                pmids = []
                assert isinstance(ea_value, list)
                # ndex:citations are in the form 'pmid:xxxxx'
                for cit in ea_value:
                    if cit.upper().startswith('PMID:'):
                        pmid = cit[5:]
                        if pmid: # Check for empty PMID strings!
                            pmids.append(pmid)
                    else:
                        logger.info("Unexpected PMID format: %s" % cit)
                ea_info['pmids'] += pmids

    def _initialize_node_attributes(self):
        self._node_attributes = _get_dict_from_list('nodeAttributes', self.cx)

    def get_agents(self):
        """Get list of grounded nodes in the network as Agents.

        Returns
        -------
        list of Agents
            Only nodes containing sufficient information to be grounded will
            be contained in this list.
        """
        return [ag for ag in self._node_agents.values()]

    def get_node_names(self):
        """Get list of all nodes in the network by name."""
        return [name for name in self._node_names.values()]

    def get_pmids(self):
        """Get list of all PMIDs associated with edges in the network."""
        pmids = []
        for ea in self._edge_attributes.values():
            edge_pmids = ea.get('pmids')
            if edge_pmids:
                pmids += edge_pmids
        return list(set(pmids))

    def get_statements(self):
        """Convert network edges into Statements.

        Returns
        -------
        list of Statements
            Converted INDRA Statements.
        """
        edges = _get_dict_from_list('edges', self.cx)
        for edge in edges:
            edge_type = edge.get('i')
            if not edge_type:
                continue
            stmt_type = _stmt_map.get(edge_type)
            if stmt_type:
                id = edge['@id']
                source_agent = self._node_agents.get(edge['s'])
                target_agent = self._node_agents.get(edge['t'])
                if not source_agent or not target_agent:
                    logger.info("Skipping edge %s->%s: %s" %
                                (self._node_names[edge['s']],
                                 self._node_names[edge['t']], edge))
                    continue
                ev = self._create_evidence(id)
                if stmt_type == Complex:
                    stmt = stmt_type([source_agent, target_agent], evidence=ev)
                else:
                    stmt = stmt_type(source_agent, target_agent, evidence=ev)
                self.statements.append(stmt)
        return self.statements

    def _create_evidence(self, edge_id):
        """Create Evidence object for a specific edge/Statement in the network.

        Parameters
        ----------
        edge_id : int
            ID of the edge in the underlying NDEx network.
        """
        pmids = None
        edge_attr = self._edge_attributes.get(edge_id)
        if edge_attr:
            pmids = edge_attr.get('pmids')
        if not pmids:
            return [Evidence(source_api='ndex',
                             source_id=self._network_info['externalId'],
                             annotations={'edge_id': edge_id})]
        else:
            evidence = []
            for pmid in pmids:
                evidence.append(
                        Evidence(source_api='ndex',
                                 source_id=self._network_info['externalId'],
                                 pmid=pmid,
                                 annotations={'edge_id': edge_id}))
            return evidence

    def get_aliases(self, node):
        cx_db_refs = {}
        node_id = node['@id']
        alias_attrs = [attr for attr in self._node_attributes if
                       attr.get('po') == node_id and
                       attr.get('n') == 'alias']
        if not alias_attrs:
            return {}
        if len(alias_attrs) > 1:
            logger.warning('More than one alias attribute for node %d' %
                           node_id)
        aliases = alias_attrs[0].get('v')
        for alias in aliases:
            if ':' not in alias:
                continue
            db_name, db_id = alias.split(':')
            db_name_mapped = cx_indra_db_map.get(db_name)
            if not db_name_mapped:
                logger.warning('DB name %s is not mapped to INDRA.' % db_name)
                continue
            cx_db_refs[db_name_mapped] = db_id
        return cx_db_refs

cx_indra_db_map = {
        'UniProt': 'UP',
        'uniprot knowledgebase': 'UP'
        }
