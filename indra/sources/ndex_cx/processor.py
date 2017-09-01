from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import objectpath
from indra.databases import uniprot_client, chebi_client, hgnc_client
from indra.literature import id_lookup
from indra.statements import *

logger = logging.getLogger('ndex_cx_processor')


_stmt_map = {
    'phosphorylates': Phosphorylation,
    'in-complex-with': Complex,
    'NEGATIVE_INFLUENCE': Inhibition,
    'POSITIVE_INFLUENCE': Activation
}


class NdexCxProcessor(object):
    def __init__(self, cx):
        self.cx = cx
        self.statements = []
        # Initialize the dict mapping node IDs to gene names
        self._node_names = {}
        self._node_agents = {}
        self._initialize_node_agents()

    def _initialize_node_agents(self):
        nodes = self.cx[5]['nodes']
        invalid_genes = []
        for node in nodes:
            id = node['@id']
            node_name = node['n']
            self._node_names[id] = node_name
            hgnc_id = hgnc_client.get_hgnc_id(node_name)
            if not hgnc_id:
                invalid_genes.append(node_name)
            else:
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                assert up_id
                self._node_agents[id] = Agent(node_name,
                                              db_refs={'HGNC': hgnc_id,
                                                       'UP': up_id})
        logger.info('Skipped invalid gene symbols: %s' %
                    ', '.join(invalid_genes))

    def get_agents(self):
        return [ag for ag in self._node_agents.values()]

    def get_node_names(self):
        return [name for name in self._node_names.values()]

    def get_statements(self):
        edges = self.cx[6]['edges']
        for edge in edges:
            edge_type = edge['i']
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
                if stmt_type == Complex:
                    stmt = stmt_type([source_agent, target_agent])
                else:
                    stmt = stmt_type(source_agent, target_agent)
                self.statements.append(stmt)
        return self.statements

    def _get_evidence(self, card):
        pmcid = card.get('pmc_id')
        ids = id_lookup(pmcid, 'pmcid')
        pmid = ids.get('pmid')
        evidence = card.get('evidence')
        all_evidence = []
        if evidence is not None:
            for text in evidence:
                e = Evidence(self.source_api, pmid=pmid, text=text)
                all_evidence.append(e)
        return all_evidence

