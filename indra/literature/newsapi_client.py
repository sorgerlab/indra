"""This module provides a client for the NewsAPI web service
(https://newsapi.org/). The web service requires an API key which is
available after registering at https://newsapi.org/account. This key can be
set as NEWSAPI_API_KEY in the INDRA config file or as an environmental
variable with the same name.

NewsAPI also requires attribution e.g. "powered by NewsAPI.org" for derived
uses.
"""
import math
import logging
import requests
from indra import has_config, get_config


logger = logging.getLogger('newsapi_client')


api_key = None
if not has_config('NEWSAPI_API_KEY'):
    logger.error('NewsAPI API key could not be found in config file or ' + \
                 'environment variable.')
else:
    api_key = get_config('NEWSAPI_API_KEY')

newsapi_url = 'https://newsapi.org/v2'


def send_request(endpoint, **kwargs):
    """Return the response to a query as JSON from the NewsAPI web service.

    The basic API is limited to 100 results which is chosen unless explicitly
    given as an argument. Beyond that, paging is supported through the "page"
    argument, if needed.

    Parameters
    ----------
    endpoint : str
        Endpoint to query, e.g. "everything" or "top-headlines"

    kwargs : dict
        A list of keyword arguments passed as parameters with the query.
        The basic ones are "q" which is the search query, "from" is a start
        date formatted as for instance 2018-06-10 and "to" is an end date
        with the same format.

    Returns
    -------
    res_json : dict
        The response from the web service as a JSON dict.
    """
    if api_key is None:
        logger.error('NewsAPI cannot be used without an API key')
        return None
    url = '%s/%s' % (newsapi_url, endpoint)
    if 'apiKey' not in kwargs:
        kwargs['apiKey'] = api_key
    if 'pageSize' not in kwargs:
        kwargs['pageSize'] = 100
    res = requests.get(url, params=kwargs)
    res.raise_for_status()
    res_json = res.json() 
    return res_json
