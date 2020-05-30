"""Client to access the Entrez Taxonomy web service."""
import requests

base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'


def _send_search_request(term):
    params = {
        'db': 'taxonomy',
        'term': term,
        'retmode': 'json'
        }
    res = requests.get(base_url, params=params)
    if not res.status_code == 200:
        return None
    return res.json().get('esearchresult')


def get_taxonomy_id(name):
    """Return the taxonomy ID corresponding to a taxonomy name.

    Parameters
    ----------
    name : str
        The name of the taxonomy entry.
        Example: "Severe acute respiratory syndrome coronavirus 2"

    Returns
    -------
    str or None
        The taxonomy ID corresponding to the given name or None
        if not available.
    """
    res = _send_search_request(name)
    idlist = res.get('idlist')
    if not idlist:
        return None
    return idlist[0]
