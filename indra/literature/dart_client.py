import json
import logging
import requests
from indra.config import get_config


logger = logging.getLogger(__name__)


dart_uname = get_config('DART_WM_USERNAME')
dart_pwd = get_config('DART_WM_PASSWORD')


dart_url = 'https://indra-ingest-pipeline-rest-1.prod.dart.worldmodelers.com' \
           '/dart/api/v1/readers/query'


def check_timestamp_dict(timestamp):
    """

    Parameters
    ----------
    timestamp : dict

    Returns
    -------
    dict
    """
    ek = {'on', 'before', 'after'}
    if sum(k in ek for k in timestamp) > 0:
        if 'on' in timestamp:
            logger.warning('Ignoring any other keys than "on"')
            ts = {'on': timestamp['on']}
        else:
            ts = {k: v for k, v in timestamp.items() if k in ek}
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
    if readers:
        pd['readers'] = readers
    if versions:
        pd['versions'] = versions
    if document_ids:
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
