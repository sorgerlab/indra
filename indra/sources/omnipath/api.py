import logging
import requests
from .processor import OmniPathProcessor

logger = logging.getLogger(__name__)


op_url = 'http://omnipathdb.org'


def process_from_web():
    """Query the omnipath web API and return an OmniPathProcessor"""
    ptm_json = _get_modifications()
    ligrec_json = _get_interactions()
    return OmniPathProcessor(ptm_json=ptm_json, ligrec_json=ligrec_json)


def _get_modifications():
    """Get all PTMs from Omnipath in JSON format.

    Returns
    -------
    JSON content for PTMs.
    """
    params = {'format': 'json',
              'fields': ['curation_effort', 'isoforms', 'references',
                         'resources', 'sources']}
    ptm_url = '%s/ptms' % op_url
    res = requests.get(ptm_url, params=params)
    if not res.status_code == 200 or not res.text:
        return None
    else:
        return res.json()


def _get_interactions(datasets=None):
    """Wrapper for calling the omnipath interactions API

    See full list of query options here:
    https://omnipathdb.org/queries/interactions

    Parameters
    ----------
    datasets
        A list of dataset names. Options are:
            dorothea, kinaseextra, ligrecextra, lncrna_mrna, mirnatarget,
            omnipath, pathwayextra, tf_mirna, tf_target, tfregulons
        Default: 'ligrecextra'

    Returns
    -------
    dict
        json of database request
    """
    interactions_url = '%s/interactions' % op_url
    params = {
        'fields': ['curation_effort', 'entity_type', 'references',
                   'resources', 'sources', 'type'],
        'format': 'json',
        'datasets': datasets or ['ligrecextra']
    }
    res = requests.get(interactions_url, params=params)
    res.raise_for_status()

    return res.json()
