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
            if inter['interaction_type'] in \
                ('adds_modification', 'removes_modification'):
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
        pass

    def _get_mod_condition(self, mod):
        mod_type = mod.get('modification_type')
        residue = mod.get('aa_code')
        position = mod.get('location')
        mc = ModCondition(mod_type, residue, position, True)
        return mc

    def _get_agent(self, participant):
        entity_type = participant.get('entity_type')
        dbid = participant.get('identifier')
        db_refs = {}
        if dbid == 'GENERIC':
            return None
        elif entity_type in ['protein', 'chemical']:
            name = participant.get('entity_text')[0]
            db_name, db_id = dbid.split(':')
            if db_name.lower() == 'uniprot':
                # TODO: get UP ID from menmonic
                db_refs['UP'] = db_id
            elif db_name.lower() == 'pubchem':
                # TODO: get ChEBI ID from PUBCHEM
                db_refs['CHEBI'] = db_id
            agent = Agent(name, db_refs=db_refs)
            return agent
        else:
            return None
        # TODO: handle other participant types

    def _get_evidence(self, card):
        pmcid = card.get('pmc_id')
        # TODO: PMCID to PMID conversion
        pmid = pmcid
        evidence = card.get('evidence')
        all_evidence = []
        if evidence is not None:
            for text in evidence:
                e = Evidence(self.source_api, pmid, text=text)
                all_evidence.append(e)
        return all_evidence

    _mod_type_map = {
            'phosphorylation': Phosphorylation,
            'dephosphorylation': Dephosphorylation,
            'ubiquitination': Ubiquitination,
            'deubiquitination': Deubiquitination
            # TODO: complete this dict
            }
