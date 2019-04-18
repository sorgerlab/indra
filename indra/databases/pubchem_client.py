import logging
import requests
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
