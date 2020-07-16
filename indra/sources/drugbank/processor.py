import logging
from xml.etree import ElementTree
from indra.statements import *
from indra.ontology.standardize import standardize_name_db_refs

logger = logging.getLogger(__name__)

drugbank_ns = {'db': 'http://www.drugbank.ca'}


class DrugbankProcessor:
    """Processor to extract INDRA Statements from DrugBank content.

    The processor assumes that an ElementTree is available which it then
    traverses to find drug-target information.

    Parameters
    ----------
    xml_tree : xml.etree.ElementTree.ElementTree
        An XML ElementTree representing DrugBank XML content.

    Attributes
    ----------
    statements : list of indra.statements.Statement
        A list of INDRA Statements that were extracted from DrugBank content.
    """
    def __init__(self, xml_tree: ElementTree.ElementTree):
        self.xml_tree = xml_tree
        self.statements = []

    def extract_statements(self):
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
            if not actions:
                actions = {'N/A'}
            for action in actions:
                stmt_type = DrugbankProcessor._get_statement_type(action)
                if not stmt_type:
                    continue
                annotations = {'drugbank_action': action}
                evs = DrugbankProcessor._get_evidences(target_element)
                for ev in evs:
                    ev.annotations = annotations
                target = DrugbankProcessor._get_target_agent(target_element)
                yield stmt_type(drug, target, evidence=evs)

    @staticmethod
    def _get_statement_type(action):
        if action in neutral_actions:
            return None
        elif action in activation_actions:
            return Activation
        elif action in inhibition_actions:
            return Inhibition
        elif action in decrease_amount_actions:
            return DecreaseAmount
        elif action in increase_amount_actions:
            return IncreaseAmount
        elif action == 'N/A':
            return Inhibition
        else:
            return None

    @staticmethod
    def _get_target_agent(target_element):
        name_tag = db_find(target_element, 'db:name')
        name = name_tag.text

        db_refs = {}

        # Get Drugbank target ID
        target_id = db_find(target_element, 'db:id').text
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
        if cas_tag is not None and cas_tag.text is not None:
            db_refs['CAS'] = cas_tag.text

        # Extract other xrefs
        for xref_tag in db_findall(drug_element, 'db:external-identifiers/'
                                   'db:external-identifier'):
            resource = db_find(xref_tag, 'db:resource').text
            identifier = db_find(xref_tag, 'db:identifier').text
            if resource == 'ChEMBL':
                db_refs['CHEMBL'] = identifier
            elif resource == 'PubChem Compound':
                db_refs['PUBCHEM'] = identifier
            elif resource == 'ChEBI':
                db_refs['CHEBI'] = identifier

        standard_name, db_refs = standardize_name_db_refs(db_refs)
        if standard_name:
            name = standard_name
        agent = Agent(name, db_refs=db_refs)
        return agent

    @staticmethod
    def _get_evidences(target_element):
        # TODO: is there a source ID we can use here?
        # TODO: is there context we can extract?
        # refs also has: textbooks, attachments
        pmids = db_findall(target_element,
                           'db:references/db:articles/db:article/db:pubmed-id')
        urls = db_findall(target_element,
                          'db:references/db:links/db:link/db:url')
        if pmids:
            evs = [Evidence(source_api='drugbank', pmid=pmid.text)
                   for pmid in pmids]
        elif urls:
            evs = [Evidence(source_api='drugbank',
                            text_refs={'URL': url.text})
                   for url in urls]
        else:
            evs = [Evidence(source_api='drugbank')]
        return evs


def db_find(element, path):
    return element.find(path, namespaces=drugbank_ns)


def db_findall(element, path):
    return element.findall(path, namespaces=drugbank_ns)


activation_actions = {'substrate', 'agonist', 'inducer', 'potentiator',
                      'stimulator', 'cofactor', 'activator', 'ligand',
                      'chaperone', 'partial agonist', 'protector',
                      'positive allosteric modulator', 'positive modulator'}

inhibition_actions = {'antagonist', 'inhibitor', 'binder', 'antibody',
                      'inactivator', 'binding', 'blocker', 'negative modulator',
                      'inverse agonist', 'neutralizer', 'weak inhibitor',
                      'suppressor', 'disruptor',
                      'inhibitory allosteric modulator'}

decrease_amount_actions = {'downregulator', 'metabolizer', 'chelator',
                           'degradation',
                           'incorporation into and destabilization'}

increase_amount_actions = {'stabilization'}

neutral_actions = {'modulator', 'other/unknown', 'unknown', 'other',
                   'regulator'}