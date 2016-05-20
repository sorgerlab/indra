import requests
import json
import time
import logging

logger = logging.getLogger('ndex')

ndex_base_url = 'http://bel2rdf.bigmech.ndexbio.org'
#ndex_base_url = 'http://52.37.175.128'

def send_request(url_suffix, params):
    """Send a request to the NDEx server.

    Parameters
    ----------
    url_suffix : str
        The API endpoint to concatenate with the base URL.
    params : dict
        A dictionary of parameters to send with the request. Parameter keys
        differ based on the type of request.

    Returns
    -------
    res_json : str
        The result of the request.
    """
    res = requests.post(ndex_base_url + url_suffix, data=json.dumps(params))
    res_json = _get_result(res)
    return res_json

def _get_result(res):
    status = res.status_code
    # If response is immediate, we get 200 
    if status == 200:
        return res.text
    # If there is a continuation of the message
    # we get status 300, handled below. 
    # Otherwise we return None.
    elif status != 300:
        return None
    task_id = res.json()['task_id']
    logger.info('NDEx task submitted...')
    time_used = 0
    try:
        while status != 200:
            res = requests.get(ndex_base_url + '/task/' + task_id)
            status = res.status_code
            if status != 200:
                time.sleep(5)
                time_used += 5
    except KeyError:
        next
        return None
    logger.info('NDEx task complete.')
    return res.text
