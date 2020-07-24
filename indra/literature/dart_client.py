import json
import logging
import requests
from datetime import datetime
from indra.config import get_config


logger = logging.getLogger(__name__)


dart_uname = get_config('DART_WM_USERNAME', failure_ok=False)
dart_pwd = get_config('DART_WM_PASSWORD', failure_ok=False)


dart_url = 'https://indra-ingest-pipeline-rest-1.prod.dart.worldmodelers.com' \
           '/dart/api/v1/readers/query'


def _check_lists(lst):
    if not isinstance(lst, (list, tuple)):
        return False
    elif any(not isinstance(s, str) for s in lst):
        logger.warning('At least one object in list is not a string')
        return False
    return True


def check_timestamp_dict(timestamp):
    """

    Parameters
    ----------
    timestamp : dict

    Returns
    -------
    dict
    """
    def _is_valid_ts(ts):
        try:
            dt = datetime.utcfromtimestamp(int(ts))
            if dt < datetime(year=1900):
                return False
        except (ValueError, OverflowError):
            return False
        return True

    ek = {'on', 'before', 'after'}
    if sum(k in ek for k in timestamp) > 0:
        if 'on' in timestamp:
            logger.warning('Ignoring any other keys than "on"')
            ts = {'on': timestamp['on']}
        else:
            ts = {k: v for k, v in timestamp.items() if k in ek and
                  _is_valid_ts(v)}
    else:
        raise ValueError(f'None of the allowed keys '
                         f'{", ".join(list(ek))} were provided')
    return ts


def jsonify_query_data(readers=None, versions=None, document_ids=None,
                       timestamp=None):
    """Check and json.dumps the metadata dictionary

    Query json structure:
    {"readers": [
        "MyAwesomeTool",
        "SomeOtherAwesomeTool"
      ],
      "versions": [
        "3.1.4",
        "1.3.3.7"
      ],
      "document_ids": [
        "qwerty1234",
        "poiuyt0987"
      ],
      "timestamp": {
        "before": {},
        "after": {},
        "on": {}}}

    Parameters
    ----------
    readers : list
    versions : list
    document_ids : list
    timestamp : dict("on"|"before"|"after",str)

    Returns
    -------
    str
        The json.dumps representation of the query metadata
    """
    if all(v is None for v in [readers, versions, document_ids, timestamp]):
        logger.warning('No query parameters were filled out')
        return ''
    pd = {}
    if readers and _check_lists(readers):
        pd['readers'] = readers
    if versions and _check_lists(versions):
        pd['versions'] = versions
    if document_ids and _check_lists(document_ids):
        pd['document_ids'] = document_ids
    if isinstance(timestamp, dict):
        pd['timestamp'] = check_timestamp_dict(timestamp)
    elif timestamp is not None:
        raise ValueError('Argument "timestamp" must be of type dict')

    return json.dumps(pd)


def query_dart_api(readers=None, versions=None, document_ids=None,
                   timestamp=None):
    """

    Parameters
    ----------
    readers : list
    versions : list
    document_ids : list
    timestamp : dict("on"|"before"|"after",str)

    Returns
    -------
    dict
        The JSON payload of the response from the DART API
    """
    query_data = jsonify_query_data(readers, versions, document_ids, timestamp)
    if not query_data:
        return {}
    res = requests.post(dart_url, data={'metadata': query_data},
                        auth=(dart_uname, dart_pwd))
    if res.status_code != 200:
        logger.warning(f'Dart Notifications Endpoint returned with status '
                       f'{res.status_code}: {res.text}')

        return {}
    return res.json()
