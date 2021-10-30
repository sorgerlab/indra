"""A client to the Biolookup web service available at http://biolookup.io/."""

from typing import Dict
import requests

URL = 'http://biolookup.io/api/lookup/'

def lookup_curie(curie: str) -> Dict:
    """Look up a CURIE in the Biolookup web service.

    Parameters
    ----------
    curie :
        The CURIE to look up.

    Returns
    -------
    :
        A dictionary containing the results of the lookup.
    """
    url = URL + curie
    response = requests.get(url)
    response.raise_for_status()
    return response.json()


def lookup(db_ns: str, db_id: str) -> dict:
    """Look up a namespace and corresponding ID in the Biolookup web service.

    Parameters
    ----------
    db_ns :
        The database namespace.
    db_id :
        The database ID.

    Returns
    -------
    :
        A dictionary containing the results of the lookup.
    """
    curie = db_ns + ':' + db_id
    return lookup_curie(curie)


def get_name(db_ns: str, db_id: str) -> Dict:
    """Return the name of a namespace and corresponding ID in the Biolookup web
    service.

    Parameters
    ----------
    db_ns :
        The database namespace.
    db_id :
        The database ID.

    Returns
    -------
    :
        The name of the entry.
    """
    res = lookup(db_ns, db_id)
    return res.get('name')
