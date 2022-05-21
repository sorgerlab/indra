__all__ = ['SemRepXmlProcessor', 'default_predicate_mappings']

from indra.ontology.standardize import get_standard_agent, \
    standardize_agent_name
from indra.statements import *


class SemRepXmlProcessor:
    def __init__(self, tree, use_gilda_grounding=False,
                 predicate_mappings=None):
        self.tree = tree
        self.entities = self.make_entity_lookup()
        self.statements = []
        self.use_gilda_grounding = use_gilda_grounding
        self.predicate_mappings = predicate_mappings if predicate_mappings \
            else default_predicate_mappings

    def process_statements(self):
        for doc in self.tree.findall('Document'):
            for utterance in doc.findall('Utterance'):
                for predication in utterance.findall('Predication'):
                    stmt = self.extract_predication(predication, utterance)
                    if stmt:
                        self.statements.append(stmt)

    def extract_predication(self, predication, utterance):
        subj = predication.find('Subject')
        obj = predication.find('Object')
        pred = predication.find('Predicate')
        negated = predication.attrib.get('negated')
        predicate_type = pred.attrib.get('type')
        stmt_type = self.predicate_mappings.get(predicate_type)
        if not stmt_type:
            return
        subj_agent = self.get_agent_from_entity(
            self.entities[subj.attrib['entityID']])
        obj_agent = self.get_agent_from_entity(
            self.entities[obj.attrib['entityID']])
        if negated:
            epi = {'negated': True}
        else:
            epi = {}
        evidence = Evidence(text=utterance.get('text'),
                            annotations={'semrep_predicate': predicate_type},
                            epistemics=epi)
        stmt = stmt_type(subj_agent, obj_agent,
                         evidence=evidence)
        return stmt

    def make_entity_lookup(self):
        return {
            entity.attrib['id']: entity
            for entity in self.tree.findall('Document/Utterance/Entity')
        }

    def get_agent_from_entity(self, entity):
        # Note: entities can be negated ("negated") and have a semantic type
        # (semtype) and score (score)
        # <Entity id="Dtest.txt.E8" cui="C3192263" name="Vemurafenib"
        # semtypes="orch,phsu" text="vemurafenib" score="851" negated="false"
        # begin="147" end="158" />
        name = entity.attrib['name']
        db_refs = {'TEXT': entity.attrib['text'],
                   'UMLS': entity.attrib['cui']}
        agent = get_standard_agent(name, db_refs)
        # We optionally add groundings from Gilda if standardization didn't
        # yield and additional references beyond UMLS.
        if self.use_gilda_grounding and set(db_refs) == {'TEXT', 'UMLS'}:
            import gilda
            matches = gilda.ground(name)
            if matches:
                db_refs[matches[0].term.db] = matches[0].term.id
                standardize_agent_name(agent, standardize_refs=True)
        return agent


# Default mappings from SemRep predicates to INDRA Statement types
default_predicate_mappings = {
    'TREATS': Inhibition,
    'INHIBITS': Inhibition,
    'PREVENTS': Inhibition,
    'DISRUPTS': Inhibition,
    'STIMULATES': Activation,
    'COMPLICATES': Activation,
    'AUGMENTS': Activation,
    'CAUSES': Activation,
    'PRODUCES': IncreaseAmount,
}