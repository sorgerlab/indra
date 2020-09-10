from indra.statements import *


class BioFactoidProcessor:
    def __init__(self, biofactoid_json):
        self.biofactoid_json = biofactoid_json
        self.statements = []

    def extract_statements(self):
        for document in self.biofactoid_json:
            self.statements += self.process_document(document)

    def process_document(self, document):
        stmts = []
        for interaction in self.find_interactions(document['elements']):
            if interaction['type'] == 'phosphorylation':
                sub = None
                enz = None
                stmt_type = None
                for entry in interaction['entries']:
                    element = self.find_element(document['elements'],
                                                entry['id'])
                    agent = self.agent_from_element(element)
                    if entry['group'] is None:
                        enz = agent
                    elif entry['group'] == 'positive':
                        sub = agent
                        stmt_type = Phosphorylation
                    elif entry['group'] == 'negative':
                        sub = agent
                        stmt_type = Dephosphorylation
                if sub and enz and stmt_type:
                    stmts.append(stmt_type(enz, sub))
        return stmts

    def agent_from_element(self, element):
        name = element['name']
        db_refs = {ref['db']: ref['id']
                   for ref in element['association']['dbXrefs']}
        return Agent(name, db_refs=db_refs)

    def find_element(self, elements, element_id):
        for element in elements:
            if element['id'] == element_id:
                return element

    def find_interactions(self, elements):
        return [e for e in elements if e['type'] in interaction_types]


# TODO: there's probably more but this is what is visible so far
interaction_types = {'binding', 'interaction', 'modification',
                     'phosphorylation', 'transcription-translation'}
