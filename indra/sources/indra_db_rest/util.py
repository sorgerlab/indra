import json
import logging
import requests

from indra import get_config
from indra.sources.indra_db_rest.exceptions import IndraDBRestAPIError

logger = logging.getLogger(__name__)


def submit_query_request(end_point, *args, **kwargs):
    """Low level function to format the query string."""
    ev_limit = kwargs.pop('ev_limit', 10)
    best_first = kwargs.pop('best_first', True)
    tries = kwargs.pop('tries', 2)
    timeout = kwargs.pop('timeout', None)
    # This isn't handled by requests because of the multiple identical agent
    # keys, e.g. {'agent': 'MEK', 'agent': 'ERK'} which is not supported in
    # python, but is allowed and necessary in these query strings.
    # TODO because we use the API Gateway, this feature is not longer needed.
    # We should just use the requests parameters dict.
    query_str = '?' + '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()
                                if v is not None]
                               + list(args))
    return submit_statement_request('get', end_point, query_str,
                                    ev_limit=ev_limit, best_first=best_first,
                                    tries=tries, timeout=timeout)


def submit_statement_request(meth, end_point, query_str='', data=None,
                             tries=2, timeout=None, **params):
    """Even lower level function to make the request."""
    full_end_point = 'statements/' + end_point.lstrip('/')
    return make_db_rest_request(meth, full_end_point, query_str, data,
                                params, tries, timeout)


def make_db_rest_request(meth, end_point, query_str, data=None, params=None,
                         tries=2, timeout=None):
    if params is None:
        params = {}

    if end_point is None:
        logger.error("Exception in submit request with args: %s"
                     % str([meth, end_point, query_str, data, params, tries]))
        raise ValueError("end_point cannot be None.")
    url_path = get_url_base(end_point)
    api_key = get_config('INDRA_DB_REST_API_KEY', failure_ok=True)
    url_path += query_str
    headers = {}
    if data:
        # This is an assumption which applies to our use cases for now, but may
        # not generalize.
        headers['content-type'] = 'application/json'
        json_data = json.dumps(data)
    else:
        json_data = None
    params['api_key'] = api_key
    logger.info('query: %s', url_path.replace(str(api_key), '[api-key]'))
    logger.info('params: %s', str(params).replace(str(api_key), '[api-key]'))
    logger.debug('headers: %s', str(headers).replace(str(api_key),
                                                     '[api-key]'))
    logger.debug('data: %s', str(data).replace(str(api_key), '[api-key]'))
    method_func = getattr(requests, meth.lower())
    while tries > 0:
        tries -= 1
        resp = method_func(url_path, headers=headers, data=json_data,
                           params=params, timeout=timeout)
        if resp.status_code == 200:
            return resp
        elif resp.status_code == 504 and tries > 0:
            logger.warning("Endpoint timed out. Trying again...")
        else:
            raise IndraDBRestAPIError(resp)


def get_url_base(end_point):
    url = get_config('INDRA_DB_REST_URL', failure_ok=False)
    url_path = url.rstrip('/') + '/' + end_point.lstrip('/')
    return url_path
