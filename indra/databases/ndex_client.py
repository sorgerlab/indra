from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import requests
import json
import time
import logging

logger = logging.getLogger('ndex')

ndex_base_url = 'http://52.37.175.128'

def send_request(ndex_service_url, params, is_json=True, use_get=False):
    """Send a request to the NDEx server.

    Parameters
    ----------
    ndex_service_url : str
        The URL of the service to use for the request.
    params : dict
        A dictionary of parameters to send with the request. Parameter keys
        differ based on the type of request.
    is_json : bool
        True if the response is in json format, otherwise it is assumed to be
        text. Default: False
    use_get : bool
        True if the request needs to use GET instead of POST.

    Returns
    -------
    res : str
        Depending on the type of service and the is_json parameter, this
        function either returns a text string or a json dict.
    """
    if use_get:
        res = requests.get(ndex_service_url, json=params)
    else:
        res = requests.post(ndex_service_url, json=params)
    status = res.status_code
    # If response is immediate, we get 200
    if status == 200:
        if is_json:
            return res.json()
        else:
            return res.text
    # If there is a continuation of the message we get status 300, handled below.
    # Otherwise we return None.
    elif status != 300:
        return None
    # In case the response is not immediate, a task ID can be used to get
    # the result.
    task_id = res.json().get('task_id')
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
    if is_json:
        return res.json()
    else:
        return res.text
