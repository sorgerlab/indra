from copy import deepcopy
from collections import defaultdict
from indra.statements import *
from indra.util import flatten


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
            groups = defaultdict(list)
            for entry in interaction['entries']:
                groups[entry['group']].append(
                    self.agent_from_element(
                        self.find_element(document['elements'], entry['id'])
                    )
                )
            groups = dict(groups)
            if interaction['type'] in mod_types:
                stmt_type = interaction_types[interaction['type']]
                enz = groups[None][0] if None in groups else None
                sub = None
                act_inh_stmt = None
                for polarity in {'positive', 'unsigned', 'negative'}:
                    if polarity in groups:
                        if polarity == 'positive':
                            act_inh_stmt = Activation
                        elif polarity == 'negative':
                            act_inh_stmt = Inhibition
                        sub = groups[polarity][0]
                        break
                if sub and enz:
                    stmt = stmt_type(deepcopy(enz),
                                     deepcopy(sub),
                                     evidence=deepcopy(ev))
                    stmts.append(stmt)
                    if act_inh_stmt:
                        stmt = act_inh_stmt(deepcopy(enz),
                                            deepcopy(sub),
                                            evidence=deepcopy(ev))
                        stmts.append(stmt)
            elif interaction['type'] == 'transcription-translation':
                subj = groups[None][0] if None in groups else None
                obj = None
                for polarity in {'positive', 'unsigned', 'negative'}:
                    if polarity in groups:
                        obj = groups[polarity][0]
                        break
                else:
                    polarity = None
                if polarity == 'positive':
                    stmt_type = IncreaseAmount
                elif polarity == 'negative':
                    stmt_type = DecreaseAmount
                if subj and obj and polarity:
                    stmt = stmt_type(deepcopy(subj),
                                     deepcopy(obj),
                                     evidence=deepcopy(ev))
                    stmts.append(stmt)
            elif interaction['type'] == 'binding':
                members = flatten(list(groups.values()))
                stmt = Complex(deepcopy(members),
                               evidence=deepcopy(ev))
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
        for ref in element['association'].get('dbXrefs', []):
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
    elif db_ns == 'RGD':
        return 'RGD', db_id
    elif db_ns == 'ChEBI':
        return 'CHEBI', 'CHEBI:%s' % db_id
    return None, None


# TODO: there's probably more but this is what is visible so far
interaction_types = {
    'binding': Complex,
    'interaction': Complex,
    'modification': AddModification,
    'phosphorylation': Phosphorylation,
    'methylation': Methylation,
    'demethylation': Demethylation,
    'ubiquitination': Ubiquitination,
    'deubiquitination': Deubiquitination,
    'transcription-translation': IncreaseAmount,
}

mod_types = {'phosphorylation', 'dephosphorylation',
             'ubiquitination', 'deubiquitination',
             'modification', 'methylation', 'demethylation'}