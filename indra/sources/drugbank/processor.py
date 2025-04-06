import logging
from xml.etree import ElementTree

from indra.databases.identifiers import ensure_chebi_prefix,\
    ensure_chembl_prefix
from indra.ontology.standardize import get_standard_agent
from indra.statements import Activation, Complex, DecreaseAmount, Evidence,\
    IncreaseAmount, Inhibition
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
                # See https://dev.drugbank.com/guides/terms/pharmacological-action
                pharm_action = db_find(target_element, 'db:known-action')
                # Skip if it's not known that it's a direct interaction
                if pharm_action.text in {'no', 'unknown'}:
                    actions = set()
                # Otherwise use the N/A action which ultimately maps to Complex
                else:
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
        elif action == 'N/A':
            return _complex
        elif action in neutral_actions:
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


activation_actions = {
    'inducer',
    'potentiator',
    'stimulator',
    'cofactor',
    'activator',
    'protector',
    'positive allosteric modulator',
    'positive modulator',
    # All agonists activate receptors, The only differences are potency,
    # how efficiently they bind and how long they stay at the receptor site
    'agonist',
    'partial agonist',
}

inhibition_actions = {
    'inhibitor',
    'antibody',
    'inactivator',
    'blocker',
    'negative modulator',
    'neutralizer',
    'weak inhibitor',
    'suppressor',
    'disruptor',
    'chelator',
    'inhibitory allosteric modulator',
    'translocation inhibitor',
    'inhibits downstream inflammation cascades',
    'nucleotide exchange blocker',
    # Antagonists can either bind to the receptor and do nothing and prevent
    # physiologic agonists to bind (which can be overcome with higher agonist
    # dosage [except in the case of irreversible antagonism which obviously
    # can't be competed with]) or can be a noncompetitive antagonist that will
    # change the structure of the active site and prevent agonist binding.
    'antagonist',
    'partial antagonist',
    # Inverse agonists act exactly the same as competitive antagonists
    # unless there is a basal physiological agonism. In which case the
    # inverse agonist will have more of an opposite effect than just a
    # pure antagonist would have.
    'inverse agonist',
}

decrease_amount_actions = {
    'downregulator',
    'metabolizer',
    'degradation',
    'incorporation into and destabilization',
    'cleavage',
    'inhibition of synthesis',
    'antisense oligonucleotide',
}

increase_amount_actions = {
    'stabilization',
    'chaperone',
    'gene replacement',
}

neutral_actions = {
    'binder',
    'binding',
    'modulator',
    'regulator',
    'substrate',
    'ligand',
    # e.g., Doxorubicin intercalates DNA to prevent transcription
    'intercalation',
    # e.g., inhibits process on a protein's aggregation (like APP or LRRK)
    'aggregation inhibitor',
    'adduct',
    'component of',
    'product of',
    'reducer',
    'oxidizer',
    # map to Ac INDRA statement?, but I'm not convinced by the idea of
    # splitting up actions
    'acetylation',
    'allosteric modulator',
    'deoxidizer',
    'cross-linking/alkylation',  # e.g. Busulfan (DB01008) alkalytes DNA
    'multitarget',
}

skip_actions = {
    'other/unknown',
    'unknown',
    'other',
}
