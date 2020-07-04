import logging
from xml.etree import ElementTree
from indra.statements import *
from indra.ontology.standardize import standardize_name_db_refs

logger = logging.getLogger(__name__)

drugbank_ns = {'db': 'http://www.drugbank.ca'}


class DrugbankProcessor:
    def __init__(self, xml_tree: ElementTree.ElementTree):
        self.xml_tree = xml_tree
        self.statements = []

    def extract_inhibitions(self):
        root = self.xml_tree.getroot()
        for drug in db_findall(root, 'db:drug'):
            for stmt in self._extract_statements_for_drug(drug):
                self.statements.append(stmt)

    @staticmethod
    def _extract_statements_for_drug(drug_element):
        drug = DrugbankProcessor._get_drug_agent(drug_element)
        for target_element in db_findall(drug_element, 'db:targets/db:target'):
            actions = {a.text for a in db_findall(target_element,
                                                  'db:actions/db:action')}
            stmt_type = DrugbankProcessor._get_statement_type(actions)
            if not stmt_type:
                continue
            annotations = {'drugbank_actions': actions}
            evs = DrugbankProcessor._get_evidences(target_element)
            for ev in evs:
                ev.annotations = annotations
            target = DrugbankProcessor._get_target_agent(target_element)
            yield stmt_type(drug, target, evidence=evs)

    @staticmethod
    def _get_statement_type(actions):
        if not actions or actions <= neutral_actions:
            return None
        if actions - neutral_actions <= inhibition_actions:
            return Inhibition
        elif actions - neutral_actions <= activation_actions:
            return Activation
        else:
            logger.error('Unhandled actions: %s' % str(actions))
            return None

    @staticmethod
    def _get_target_agent(target_element):
        name_tag = db_find(target_element, 'db:name')
        name = name_tag.text
        assert name is not None

        db_refs = {}

        # Get Drugbank target ID
        target_id = db_find(target_element, 'db:id').text
        assert target_id
        db_refs['DRUGBANKV4.TARGET'] = target_id

        # Extract other xrefs
        for xref_tag in db_findall(target_element, 'db:polypeptide/'
                                   'db:external-identifiers/'
                                   'db:external-identifier'):
            resource = db_find(xref_tag, 'db:resource').text
            identifier = db_find(xref_tag, 'db:identifier').text
            if resource == 'HUGO Gene Nomenclature Committee (HGNC)':
                db_refs['HGNC'] = identifier[5:]
            elif resource == 'UniProtKB':
                db_refs['UP'] = identifier
        standard_name, db_refs = standardize_name_db_refs(db_refs)
        if standard_name:
            name = standard_name
        agent = Agent(name, db_refs=db_refs)
        return agent

    @staticmethod
    def _get_drug_agent(drug_element):
        name_tag = db_find(drug_element, 'db:name')
        name = name_tag.text
        assert name is not None

        db_refs = {}

        # Extract the DrugBank ID
        drugbank_id_tags = db_findall(drug_element, 'db:drugbank-id')
        # We do a sort here because sometimes there's more than one
        # DrugBank ID and we choose the "smaller" one here
        drugbank_id = sorted([di.text for di in drugbank_id_tags
                              if di.text.startswith('DB')])[0]
        db_refs['DRUGBANK'] = drugbank_id

        # Extract CAS ID
        cas_tag = db_find(drug_element, 'db:cas-number')
        if cas_tag is not None:
            db_refs['CAS'] = cas_tag.text

        # Extract other xrefs
        for xref_tag in db_findall(drug_element, 'db:external-identifiers/'
                                   'db:external-identifier'):
            resource = db_find(xref_tag, 'db:resource').text
            identifier = db_find(xref_tag, 'db:identifier').text

        standard_name, db_refs = standardize_name_db_refs(db_refs)
        if standard_name:
            name = standard_name
        agent = Agent(name, db_refs=db_refs)
        return agent

    @staticmethod
    def _get_evidences(target_element):
        # TODO: is there a source ID we can use here?
        # TODO: is there context we can extract?
        # refs also has: textbooks, links, attachments
        pmids = db_findall(target_element,
                           'db:references/db:articles/db:article/db:pubmed-id')
        evs = [Evidence(source_api='drugbank', pmid=pmid.text)
               for pmid in pmids]
        return evs


def db_find(element, path):
    return element.find(path, namespaces=drugbank_ns)


def db_findall(element, path):
    return element.findall(path, namespaces=drugbank_ns)


activation_actions = {'substrate', 'agonist', 'inducer', 'potentiator',
                      'stimulator', 'cofactor', 'activator', 'ligand',
                      'chaperone', 'partial agonist', 'protector',
                      'positive allosteric modulator'}

inhibition_actions = {'antagonist', 'inhibitor', 'binder', 'antibody',
                      'inactivator', 'binding', 'blocker', 'negative modulator',
                      'inverse agonist', 'neutralizer', 'weak inhibitor',
                      'suppressor', 'disruptor'}

decrease_amount_actions = {'downregulator'}


neutral_actions = {'modulator', 'other/unknown', 'unknown', 'other',
                   'regulator'}