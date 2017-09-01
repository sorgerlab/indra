from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import objectpath
from indra.databases import uniprot_client, chebi_client, hgnc_client
from indra.literature import id_lookup
from indra.statements import *

logger = logging.getLogger('ndex_cx_processor')

class NdexCxProcessor(object):
    def __init__(self, cx):
        self.cx = cx
        self.statements = []
        # Initialize the dict mapping node IDs to gene names
        self._node_map = {}
        self._initialize_node_agents()

    def _initialize_node_agents(self):
        nodes = self.cx[5]['nodes']
        invalid_genes = []
        for node in nodes:
            id = node['@id']
            node_name = node['n']
            hgnc_id = hgnc_client.get_hgnc_id(node_name)
            if not hgnc_id:
                invalid_genes.append(node_name)
            else:
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                assert up_id
                self._node_map[id] = Agent(node_name,db_refs={'HGNC': hgnc_id,
                                                              'UP': up_id})
        logger.info('Skipped invalid gene symbols: %s' %
                    ', '.join(invalid_genes))

    def get_statements(self):
        edges = self.cx[6]['edges']
        for edge in edges:
            id = edge['@id']
            edge_type = edge['i']
            source_node = edge['s']
            target_node = edge['t']

    def _get_modifications(self): 

        for card in self.index_cards:
            inter = card.get('interaction')
            if inter['interaction_type'] not in \
                ('adds_modification', 'removes_modification'):
                continue
            enz = self._get_agent(inter.get('participant_a'))
            sub = self._get_agent(inter.get('participant_b'))
            if sub is None:
                continue
            mods = inter.get('modifications')
            mcs = [self._get_mod_condition(mod) for mod in mods]
            ev = self._get_evidence(card)
            for mc in mcs:
                stmt_class = modtype_to_modclass.get(mc.mod_type)
                if stmt_class is None:
                    print('%s not found in mod type map' % mc.mod_type)
                    continue
                stmt = stmt_class(enz, sub, mc.residue, mc.position,
                                  evidence=ev)
                self.statements.append(stmt)

    def get_binds(self): 
        for card in self.index_cards:
            inter = card.get('interaction')
            if inter['interaction_type'] != 'binds':
                continue
            a1 = self._get_agent(inter.get('participant_a'))
            a2 = self._get_agent(inter.get('participant_b'))
            if a1 is None or a2 is None:
                continue
            ev = self._get_evidence(card)
            stmt = Complex([a1, a2], evidence=ev)
            self.statements.append(stmt)

    def get_complexes(self):
        for card in self.index_cards:
            inter = card.get('interaction')
            if inter['interaction_type'] != 'complexes_with':
                continue
            ev = self._get_evidence(card)
            participant = inter.get('participant_a')
            entities = participant.get('entities')
            members = []
            for entity in entities:
                agent = self._get_agent(entity)
                if agent is not None:
                    members.append(agent)
            if len(members) < 2:
                continue
            stmt = Complex(members, evidence=ev)
            self.statements.append(stmt)

    def get_increase_decrease(self):
        for card in self.index_cards:
            inter = card.get('interaction')
            if inter['interaction_type'] not in ('decreases', 'increases'):
                continue
            ev = self._get_evidence(card)
            participant = inter.get('participant_a')
            controller = self._get_agent(participant)
            process = inter.get('participant_b')

    def get_translocates(self):
        for card in self.index_cards:
            inter = card.get('interaction')
            if inter['interaction_type'] != 'translocates':
                continue
            ev = self._get_evidence(card)
            participant = inter.get('participant_b')
            agent = self._get_agent(participant)
            from_location = inter.get('from_location_id')
            to_location = inter.get('to_location_id')
            stmt = Translocation(agent, from_location, to_location,
                                 evidence=ev)
            self.statements.append(stmt)

    def _get_agent(self, participant):
        dbid = participant.get('identifier')
        text = participant.get('entity_text')[0]

        if dbid == 'GENERIC':
            if not text:
                return None
            else:
                return Agent(text)

        db_refs = {}
        entity_type = participant.get('entity_type')
        if entity_type in ['protein', 'chemical', 'gene']:
            # TODO: standardize name here
            name = participant.get('entity_text')[0]
            db_refs['TEXT'] = text
            if dbid:
                db_name, db_id = dbid.split(':')
                if db_name.lower() == 'uniprot':
                    uniprot_id = uniprot_client.get_id_from_mnemonic(db_id)
                    db_refs['UP'] = uniprot_id
                elif db_name.lower() == 'pubchem':
                    chebi_id = chebi_client.get_chebi_id_from_pubchem(db_id)
                    db_refs['CHEBI'] = chebi_id
                elif db_name.lower() == 'hgnc':
                    db_refs['HGNC'] = db_id
        elif entity_type == 'protein_family':
            name = text
        else:
            return None
        # TODO: handle other participant types
        agent = Agent(name, db_refs=db_refs)

        features = participant.get('features')
        if features:
            for feature in features:
                feature_type = feature.get('feature_type')
                if feature_type == 'modification_feature':
                    mc = self._get_mod_condition(feature)
                    agent.mods.append(mc)
                elif feature_type == 'binding_feature':
                    bc = self._get_bound_condition(feature)
                    agent.bound_conditions.append(bc)
                elif feature_type == 'mutation_feature':
                    mc = self._get_mut_condition(feature)
                    agent.mutations.append(mc)
                elif feature_type == 'location_feature':
                    agent.location = feature.get('location')
        not_features = participant.get('features')
        if not_features:
            for feature in not_features:
                feature_type = feature.get('feature_type')
                if feature_type == 'modification_feature':
                    mc = self._get_mod_condition(feature)
                    mc.is_modified = False
                    agent.mods.append(mc)
                elif feature_type == 'binding_feature':
                    bc = self._get_bound_condition(feature)
                    bc.is_bound = False
                    agent.bound_conditions.append(bc)
        return agent

    def _get_mod_condition(self, mod):
        mod_type = mod.get('modification_type')
        residue = mod.get('aa_code')
        position = mod.get('location')
        mc = ModCondition(mod_type, residue, position, True)
        return mc

    def _get_bound_condition(self, feature):
        bound_to = feature.get('bound_to')
        agent = self._get_agent(bound_to)
        bc = BoundCondition(agent, True)
        return bc

    def _get_mut_condition(self, mod):
        mod_type = mod.get('modification_type')
        from_residue = mod.get('from_aa')
        to_residue = mod.get('to_aa')
        position = mod.get('location')
        mc = MutCondition(position, from_residue, to_residue)
        return mc

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

