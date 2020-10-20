import logging
import requests

from indra.databases import uniprot_client, hgnc_client
from indra.statements.validate import validate_text_refs
from indra.statements import Phosphorylation, Evidence, Agent

from .phosphoelm_mapping import phosphoelm_mapping

gilda_url = 'http://grounding.indra.bio/ground'
logger = logging.getLogger(__name__)


class PhosphoElmProcessor(object):
    """Processes data dumps from the phospho.ELM database.

    See http://phospho.elm.eu.org/dataset.html

    Parameters
    ----------
    phosphoelm_data : list[dict]
        JSON compatible list of entries from a phospho.ELM data dump

    Attributes
    ----------
    statements : list[indra.statements.Phosphorylation]
        A list of the phosphorylation statements produced by the entries
        in phosphoelm_data
    """
    def __init__(self, phosphoelm_data):
        self.statements = []
        self._phosphoelm_data = phosphoelm_data

    def process_phosphorylations(self, skip_empty=True):
        """Create Phosphorylation statements from phosphoelm_data

        Parameters
        ----------
        skip_empty : bool
            Default: True. If False, also create statements when upstream
            kinases in entry['kinases'] are not known.
        """
        for entry in self._phosphoelm_data:
            if entry['species'].lower() != 'homo sapiens' or\
                    skip_empty and not entry['kinases']:
                # Skip entries without any kinases or if species is other
                # than human.
                continue
            # Entries:
            # 'acc': '<UP ID>', <-- substrate
            # 'sequence': '<protein sequence>',
            # 'position': '<sequence position>',
            # 'code': '<phosphorylated residue>',
            # 'pmids': '<pmid>',
            # 'kinases': '<responsible kinase>', <-- enzyme
            # 'source': 'HTP|LTP',
            # 'species': '<species name in latin>',
            # 'entry_date': 'yyyy-mm-dd HH:MM:SS.mmmmmm'
            substrate = _agent_from_id(entry['acc'])
            enzyme = _agent_from_str(entry['kinases'])

            # Skip if enz is None instead of an Agent (only when we skip
            # empty kinase entries)
            if skip_empty and enzyme is None:
                continue

            pmid = entry['pmids']
            if not validate_text_refs({'PMID': pmid}):
                pmid = None

            # Build evidence, add statement
            evidence = Evidence(
                source_api='phosphoelm',
                pmid=pmid,
                annotations={
                    'data_source': entry.get('source'),
                    'phosphoelm_substrate_id': entry['acc'],
                    'phosphoelm_kinase_name': entry.get('kinases'),
                    'entry_date': entry['entry_date'],
                    'sequence': entry['sequence']
                }
            )
            self.statements.append(Phosphorylation(
                enz=enzyme,
                sub=substrate,
                residue=entry['code'],
                position=entry['position'],
                evidence=evidence)
            )


def _agent_from_id(db_id):
    # There are some Ensembl protein IDs which we currently can't normalize
    # to anything else (unlike ENSG).
    if db_id.startswith('ENSP'):
        db_refs = {'ENSEMBL': db_id}
        name = db_id
    # All other entries are UniProt IDs
    else:
        name = uniprot_client.get_gene_name(db_id)
        if not name:
            return None
        if '-' in db_id:
            up_base = db_id.split('-')[0]
            db_refs = {'UP': up_base, 'UPISO': db_id}
        else:
            db_refs = {'UP': db_id}
        hgnc_id = uniprot_client.get_hgnc_id(db_id)
        if hgnc_id:
            db_refs['HGNC'] = hgnc_id
    return Agent(name, db_refs=db_refs)


def _agent_from_str(txt):
    """Return a grounded Agent from the name of an entity.

    Parameters
    ----------
    txt : str
        A string representing an entity

    Returns
    -------
    ag : indra.statements.Agent
        A grounded INDRA Agent corresponding to the provided entity text.
    """
    # Check if kinase is in the hard coded mapping
    if txt in phosphoelm_mapping:
        name = txt
        ns, _id = phosphoelm_mapping[txt]
        # If None is hard coded in the map it means we should skip this
        if ns is None:
            return None
    else:
        term = _gilda_grounder(txt)

        if term is None:
            name = txt
            ns, _id = None, None
        else:
            name = term['entry_name']
            ns, _id = term['db'], term['id']

    db_refs = {'TEXT': txt}
    if ns is not None and _id is not None:
        db_refs[ns] = _id
    if ns == 'HGNC':
        up_id = hgnc_client.get_uniprot_id(_id)
        if up_id:
            db_refs['UP'] = up_id
        name = hgnc_client.get_hgnc_name(_id)
    elif ns == 'FPLX':
        name = _id
    return Agent(name, db_refs=db_refs)


def _gilda_grounder(txt):
    # Pre-process text for grounding
    txt = txt.replace('_group', '')
    txt = txt.replace('_', '-')
    txt = txt.split('/')[0]
    res = requests.post(gilda_url, json={'text': txt})
    if res.status_code != 200:
        logger.warning('Gilda service responded with status code %d' %
                       res.status_code)
        return None
    rj = res.json()
    if not rj:
        return None
    top_term = rj[0]['term']
    return top_term

