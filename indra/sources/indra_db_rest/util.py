import json
import logging
from io import StringIO
from contextlib import contextmanager

import requests

from indra import get_config
from indra.sources.indra_db_rest.exceptions import IndraDBRestAPIError


class RecordableLogger:
    def __init__(self, name):
        self.__logger = logging.getLogger(name)

        # Set up redirection of the logs surrounding requests.
        self.__logstream = StringIO()
        self.__logstream_handler = logging.StreamHandler(self.__logstream)
        fmt = "%(levelname)s: [%(asctime)s] %(name)s - %(message)s"
        self.__logstream_handler.setFormatter(logging.Formatter(fmt))
        self.__quieted = False

    def __getattr__(self, item):
        # If the attribute could not be retrieved by standard means, try passing
        # the method on to the backend logger.
        return getattr(self.__logger, item)

    def get_quiet_logs(self):
        """Get the logstream string for the quieted request logs."""
        return self.__logstream.getvalue()

    def quiet(self):
        """Stop printing logging messages to stdout/stderr.

        The log messages are preserved, and can still be accessed using the
        `get_quiet_logs` method.
        """
        if not self.__quieted:
            self.__logger.addHandler(self.__logstream_handler)
            self.__logger.propagate = False
            self.__quieted = True

    def unquiet(self):
        """Allow the logs to print to stdout/stderr.

        Log messages will no longer be stored in the StringIO buffer.
        """
        if self.__quieted:
            self.__logger.removeHandler(self.__logstream_handler)
            self.__logger.propagate = True
            self.__quieted = False

    @contextmanager
    def quieted(self):
        self.quiet()
        try:
            yield
        finally:
            self.unquiet()


logger = RecordableLogger(__name__)


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


def make_db_rest_request(meth, end_point, query_str='', data=None, params=None,
                         tries=2, timeout=None, api_key=None):
    if params is None:
        params = {}

    if end_point is None:
        logger.error("Exception in submit request with args: %s"
                     % str([meth, end_point, query_str, data, params, tries]))
        raise ValueError("end_point cannot be None.")
    url_path = get_url_base(end_point)
    if api_key is None:
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

    def remove_api_key(s):
        if api_key:
            return s.replace(str(api_key), '[api-key]')

    logger.info(f'query: {remove_api_key(url_path)}')
    logger.info(f'params: {remove_api_key(str(params))}')
    logger.info(f'data: {remove_api_key(str(data))}')
    logger.debug(f'headers: {remove_api_key(str(headers))}')
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


def jsonify_args(d):
    new_d = d.copy()
    for key, val in d.items():
        if isinstance(val, set):
            new_d[key] = list(val)
    return new_d


def get_url_base(end_point):
    url = get_config('INDRA_DB_REST_URL', failure_ok=False)
    url_path = url.rstrip('/') + '/' + end_point.lstrip('/')
    return url_path
