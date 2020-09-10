import copy
from indra.statements import *


class BioFactoidProcessor:
    def __init__(self, biofactoid_json):
        self.biofactoid_json = biofactoid_json
        self.statements = []

    def extract_statements(self):
        for document in self.biofactoid_json:
            self.statements += self.process_document(document)

    def process_document(self, document):
        ev = self.get_evidence(document)
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
                    stmt = stmt_type(enz, sub, evidence=copy.deepcopy(ev))
                    stmts.append(stmt)
        return stmts

    def agent_from_element(self, element):
        name = element['name']
        db_refs = {}
        mapped_ns, mapped_id = \
            process_db_refs(element['association']['namespace'],
                            element['association']['id'])
        if mapped_ns and mapped_id:
            db_refs[mapped_ns] = mapped_id
        for ref in element['association']['dbXrefs']:
            mapped_ns, mapped_id = process_db_refs(ref['db'], ref['id'])
            if mapped_ns and mapped_id:
                db_refs[mapped_ns] = mapped_id
        return Agent(name, db_refs=db_refs)

    def find_element(self, elements, element_id):
        for element in elements:
            if element['id'] == element_id:
                return element

    def find_interactions(self, elements):
        return [e for e in elements if e['type'] in interaction_types]

    def get_evidence(self, document):
        text_refs = self.get_text_refs(document)
        pmid = text_refs.get('PMID')
        annotations = {
            'biofactoid_document': document['id'],
            'created_date': document.get('createdDate'),
            'lsatEditedDate': document.get('lastEditedDate'),
        }
        ev = Evidence(source_api='biofactoid', pmid=pmid, text_refs=text_refs,
                      text=document['text'], annotations=annotations)
        return ev

    def get_text_refs(self, document):
        text_refs = {}
        cit = document['citation']
        if 'pmid' in cit:
            text_refs['PMID'] = cit['pmid']
        if 'doi' in cit:
            text_refs['DOI'] = cit['doi']
        pmd = document['article']['PubmedData']
        for entry in pmd['ArticleIdList']:
            if entry['IdType'] == 'pubmed':
                text_refs['PMID'] = entry['id']
            elif entry['IdType'] == 'doi':
                text_refs['DOI'] = entry['id']
            elif entry['IdType'] == 'pii':
                text_refs['PII'] = entry['id']
            elif entry['IdType'] == 'pmc':
                text_refs['PMCID'] = entry['id']
        return text_refs


def process_db_refs(db_ns, db_id):
    if db_ns == 'HGNC':
        return 'HGNC', db_id.replace('HGNC:', '')
    elif db_ns == 'Ensembl':
        return 'ENSEMBL', db_id
    elif db_ns == 'ncbi':
        return 'EGID', db_id
    elif db_ns == 'MGI':
        return 'MGI', db_id.replace('MGI:', '')
    return None, None


# TODO: there's probably more but this is what is visible so far
interaction_types = {'binding', 'interaction', 'modification',
                     'phosphorylation', 'transcription-translation'}
