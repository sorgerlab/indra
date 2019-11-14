import logging
import requests
from .processor import OmniPathProcessor


logger = logging.getLogger("omnipath")


op_url = 'http://omnipathdb.org'


def process_from_web():
    ptm_json = _get_modifications()
    return OmniPathProcessor(ptm_json)

def _get_modifications():
    """Get all PTMs from Omnipath in JSON format.

    Returns
    -------
    JSON content for PTMs.
    """
    #params = {'format': 'json', 'fields':['sources', 'references']}
    params = {'format': 'json', 'fields':['sources']}
    ptm_url = '%s/ptms' % op_url
    res = requests.get(ptm_url, params=params)
    if not res.status_code == 200 or not res.text:
        return None
    else:
        return res.json()

