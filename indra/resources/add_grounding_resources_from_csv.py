import logging
import requests
from multiprocessing import Lock
from multiprocessing.pool import ThreadPool
from indra.databases.identifiers import get_identifiers_url


logger = logging.getLogger(__name__)

lock = Lock()


def validate_identifiers(identifiers_list, n_threads=1):
    """Batch check if db_name, db_id pairs are valid identifiers

    Parameters
    ----------
    identifiers_list : list
        List of tuples of strings of the form (db_name, db_id)

    Returns
    -------
    list
        List with same number of entries as identifiers_list. Value is True
        if the corresponding entry in identifiers_list corresponds to a valid
        identifier, otherwise False.
    """
    unique_identifiers = list(set(identifiers_list))
    urls = [get_identifiers_url(db_name, db_id)
            for db_name, db_id in unique_identifiers]
    with ThreadPool(n_threads) as pool:
        validity_list = pool.map(_validate_identifier_url, urls)
    validity_map = {(db_name, db_id): validity
                    for (db_name, db_id), validity in
                    zip(unique_identifiers, validity_list)}
    return [validity_map[(db_name, db_id)]
            for db_name, db_id in identifiers_list]


def _validate_identifier_url(url):
    if url is None:
        return False
    response = requests.get(url)
    if response.status_code == 200:
        return True
    elif response.status_code in [400, 404, 410]:
        return False
    else:
        with lock:
            logger.warning('Unable to validate url: %s due to error %s'
                           % (url, response.status_code))
        return None
