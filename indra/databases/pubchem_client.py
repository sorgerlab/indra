import logging
import requests
from typing import List
from functools import lru_cache


pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'


logger = logging.getLogger(__name__)


@lru_cache(maxsize=5000)
def get_inchi_key(pubchem_cid):
    """Return the InChIKey for a given PubChem CID.

    Parameters
    ----------
    pubchem_cid : str
        The PubChem CID whose InChIKey should be returned.

    Returns
    -------
    str
        The InChIKey corresponding to the PubChem CID.
    """
    url = '%s/compound/cid/%s/property/InChIKey/TXT' % \
        (pubchem_url, pubchem_cid)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Could not retrieve InChIKey for %s' % pubchem_cid)
        return None
    return res.text.strip()


@lru_cache(maxsize=5000)
def get_json_record(pubchem_cid):
    """Return the JSON record of a given PubChem CID.

    Parameters
    ----------
    pubchem_cid : str
        The PubChem CID whose record should be returned.

    Returns
    -------
    dict
        The record deserialized from JSON.
    """
    url = (pubchem_url + '_view/data/compound/%s/JSON/') % pubchem_cid
    res = requests.get(url)
    return res.json()


def get_preferred_compound_ids(pubchem_cid):
    """Return a list of preferred CIDs for a given PubChem CID.

    Parameters
    ----------
    pubchem_cid : str
        The PubChem CID whose preferred CIDs should be returned.

    Returns
    -------
    list of str
        The list of preferred CIDs for the given CID. If there are no
        preferred CIDs for the given CID then an empty list is returned.
    """
    record = get_json_record(pubchem_cid)
    sections = record['Record']['Section']
    pref_ids = set()
    for section in sections:
        if section['TOCHeading'] == 'Preferred Compound':
            for pref_cpd in section['Information']:
                pref_ids |= set(pref_cpd['Value']['Number'])
    pref_ids = sorted([str(pid) for pid in pref_ids])
    return pref_ids


def get_pmids(pubchem_cid: str) -> List[str]:
    """Return xref PMIDs for a given PubChem CID.

    Parameters
    ----------
    pubchem_cid :
        The PubChem CID whose xref PubMedID's will be returned.

    Returns
    -------
    list of str
        PubMedIDs corresponding to the given PubChem CID. If none present,
        an empty list is returned.
    """
    url = '%s/compound/cid/%s/xrefs/PubMedID/JSON' % \
          (pubchem_url, pubchem_cid)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Could not retrieve PMIDs for %s' % pubchem_cid)
        return []
    res_json = res.json()
    pmids_list = [str(pmid) for pmid in res_json['InformationList']['Information'][0]['PubMedID']]
    return pmids_list
