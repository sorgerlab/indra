import logging
from xml.etree import ElementTree

from indra.databases.identifiers import ensure_chebi_prefix, ensure_chembl_prefix
from indra.ontology.standardize import get_standard_agent
from indra.statements import Activation, Complex, DecreaseAmount, Evidence, IncreaseAmount, Inhibition
from indra.statements.validate import assert_valid_db_refs

logger = logging.getLogger(__name__)

drugbank_ns = {'db': 'http://www.drugbank.ca'}

_UNHANDLED = set()


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
        if action in activation_actions:
            return Activation
        elif action in inhibition_actions:
            return Inhibition
        elif action in decrease_amount_actions:
            return DecreaseAmount
        elif action in increase_amount_actions:
            return IncreaseAmount
        elif action in neutral_actions or action == 'N/A':
            return _complex
        elif action in skip_actions:
            return None
        elif action not in _UNHANDLED:
            _UNHANDLED.add(action)
            logger.warning('unhandled DrugBank action: %s', action)

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
        return get_standard_agent(name, db_refs=db_refs)

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
                db_refs['CHEMBL'] = ensure_chembl_prefix(identifier)
            elif resource == 'PubChem Compound':
                db_refs['PUBCHEM'] = identifier
            elif resource == 'ChEBI':
                db_refs['CHEBI'] = ensure_chebi_prefix(identifier)
        assert_valid_db_refs(db_refs)
        return get_standard_agent(name, db_refs)

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


def _complex(a, b, evidence):
    return Complex([a, b], evidence=evidence)


activation_actions = {'inducer', 'potentiator',
                      'stimulator', 'cofactor', 'activator',
                      'protector',
                      'positive allosteric modulator', 'positive modulator'}

inhibition_actions = {'inhibitor', 'binder', 'antibody',
                      'inactivator', 'binding', 'blocker', 'negative modulator',
                      'neutralizer', 'weak inhibitor',
                      'suppressor', 'disruptor', 'chelator',
                      'inhibitory allosteric modulator'}

decrease_amount_actions = {
    'downregulator',
    'metabolizer',
    'degradation',
    'incorporation into and destabilization',
    'cleavage',
    'inhibition of synthesis',
}

increase_amount_actions = {'stabilization', 'chaperone'}

neutral_actions = {
    'modulator',
    'regulator',
    'antagonist',
    'substrate',
    'agonist',
    'ligand',
    'intercalation',  # e.g., Doxorubicin intercalates DNA to prevent transcription
    'inverse agonist',
    'aggregation inhibitor',  # e.g., inhibits process on a protein's aggregation (like APP or LRRK)
    'partial agonist',
    'partial antagonist',
    'antisense oligonucleotide',
    'adduct',
    'component of',
    'product of',
    'reducer',
    'oxidizer',
    'acetylation',  # map to Ac INDRA statement?, but I'm not convinced by the idea of splitting up actions
}

skip_actions = {
    'other/unknown',
    'unknown',
    'other',
}
