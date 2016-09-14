import json
import objectpath
from indra.statements import *

class IndexCardProcessor(object):
    def __init__(self, index_cards, source_api):
        self.index_cards = index_cards
        self.source_api = source_api
        self.statements = []

    def get_modifications(self): 
        for card in self.index_cards:
            inter = card.get('interaction')
            if inter['interaction_type'] not in \
                ('adds_modification', 'removes_modification'):
                continue
            enz = self._get_agent(inter.get('participant_a'))
            sub = self._get_agent(inter.get('participant_b'))
            mods = inter.get('modifications')
            mcs = [self._get_mod_condition(mod) for mod in mods]
            ev = self._get_evidence(card)
            for mc in mcs:
                stmt_class = self._mod_type_map.get(mc.mod_type)
                stmt = stmt_class(enz, sub, mc.residue, mc.position,
                                  evidence=ev)
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
                members.append(agent)
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
        if dbid == 'GENERIC':
            return None

        db_refs = {}
        entity_type = participant.get('entity_type')
        if entity_type in ['protein', 'chemical']:
            # TODO: standardize name here
            name = participant.get('entity_text')[0]
            db_name, db_id = dbid.split(':')
            if db_name.lower() == 'uniprot':
                # TODO: get UP ID from menmonic
                db_refs['UP'] = db_id
            elif db_name.lower() == 'pubchem':
                # TODO: get ChEBI ID from PUBCHEM
                db_refs['CHEBI'] = db_id
            db_refs['TEXT'] = participant.get('entity_text')[0]
        elif entity_type == 'protein_family':
            name = participant.get('entity_text')[0]
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
        # TODO: PMCID to PMID conversion
        pmid = pmcid
        evidence = card.get('evidence')
        all_evidence = []
        if evidence is not None:
            for text in evidence:
                e = Evidence(self.source_api, pmid=pmid, text=text)
                all_evidence.append(e)
        return all_evidence

    _mod_type_map = {
            'phosphorylation': Phosphorylation,
            'dephosphorylation': Dephosphorylation,
            'ubiquitination': Ubiquitination,
            'deubiquitination': Deubiquitination,
            'acetylation': Acetylation,
            'deacetylation': Deacetylation,
            'farnesylation': Farnesylation,
            'glycosylation': Glycosylation,
            'deglycosylation': Deglycosylation,
            'hydorxylation': Hydroxylation,
            'dehydroxylatoni': Dehydroxylation,
            'sumoylation': Sumoylation,
            'desumoylation': Desumoylation
            }
