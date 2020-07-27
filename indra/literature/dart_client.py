"""A client for accessing reader output from the DART system."""

import json
import logging
import requests
from datetime import datetime
from indra.config import get_config


logger = logging.getLogger(__name__)


dart_uname = get_config('DART_WM_USERNAME')
dart_pwd = get_config('DART_WM_PASSWORD')

dart_base_url = 'https://indra-ingest-pipeline-rest-1.prod.dart' \
                '.worldmodelers.com/dart/api/v1/readers'
meta_endpoint = dart_base_url + '/query'
downl_endpoint = dart_base_url + '/download/'


def get_documents(readers=None, versions=None, document_ids=None,
                  timestamp=None):
    """Return a list of contents constrained by the provided parameters

    Parameters
    ----------
    readers : list
        A list of reader names
    versions : list
        A list of versions to match with the reader name(s)
    document_ids : list
        A list of document identifiers
    timestamp : dict("on"|"before"|"after",str)
        The timestamp string must of format "yyyy-mm-dd" or "yyyy-mm-dd
        hh:mm:ss" (only for "before" and "after").

    Returns
    -------
    dict(str, str)
        A dict of document content keyed by document id
    """
    metdata_json = get_reader_outputs(readers=readers, versions=versions,
                                      document_ids=document_ids,
                                      timestamp=timestamp)
    # Loop document keys and get documents
    documents = {}
    if metdata_json:
        for doc_key in metdata_json:
            doc = requests.get(url=downl_endpoint + doc_key)
            documents[doc_key] = doc
    else:
        logger.warning('Empty meta data json returned')
    return documents


def get_reader_outputs(readers=None, versions=None, document_ids=None,
                       timestamp=None):
    """Return reader output metadata by querying the DART API

    Query json structure:
        {"readers": ["MyAwesomeTool", "SomeOtherAwesomeTool"],
        "versions": ["3.1.4", "1.3.3.7"],
        "document_ids": ["qwerty1234", "poiuyt0987"],
        "timestamp": {"before": "yyyy-mm-dd"|"yyyy-mm-dd hh:mm:ss",
        "after": "yyyy-mm-dd"|"yyyy-mm-dd hh:mm:ss",
        "on": "yyyy-mm-dd"}}

    Parameters
    ----------
    readers : list
        A list of reader names
    versions : list
        A list of versions to match with the reader name(s)
    document_ids : list
        A list of document identifiers
    timestamp : dict("on"|"before"|"after",str)
        The timestamp string must of format "yyyy-mm-dd" or "yyyy-mm-dd
        hh:mm:ss" (only for "before" and "after").

    Returns
    -------
    dict
        The JSON payload of the response from the DART API
    """
    if not dart_uname:
        raise ValueError('DART_WM_USERNAME is not configured.')
    if not dart_pwd:
        raise ValueError('DART_WM_PASSWORD is not configured.')
    query_data = _jsonify_query_data(readers, versions, document_ids, timestamp)
    if not query_data:
        return {}
    res = requests.post(meta_endpoint, data={'metadata': query_data},
                        auth=(dart_uname, dart_pwd))
    res.raise_for_status()
    return res.json()


def _check_lists(lst):
    if not isinstance(lst, (list, tuple)):
        return False
    elif any(not isinstance(s, str) for s in lst):
        logger.warning('At least one object in list is not a string')
        return False
    return True


def _check_timestamp_dict(ts_dict):
    """Check the timestamp dict

    Parameters
    ----------
    ts_dict : dict
        Timestamp should be of format "yyyy-mm-dd". "yyyy-mm-dd hh:mm:ss"
        is allowed as well for the keys "before" and "after".

    Returns
    -------
    dict
    """
    def _is_valid_ts(k, tstr):
        """
        %Y - Year as Zero padded decimal
        %m - month as zero padded number
        %d - day as zero padded number
        %H - 24h hour as zero padded number
        %M - minute as zero padded number
        %S - second as zero padded number
        """
        ts_fmt = '%Y-%m-%d'
        ts_long_fmt = '%Y-%m-%d %H:%M:%S'
        if k == 'on':
            dt = datetime.strptime(tstr, ts_fmt)
        else:
            try:
                dt = datetime.strptime(tstr, ts_long_fmt)
            except ValueError:
                try:
                    dt = datetime.strptime(tstr, ts_fmt)
                except ValueError as err:
                    raise ValueError(
                        f'Timestamp "{tstr}" is not in a valid format. '
                        f'Format must be "%Y-%m-%d" or "%Y-%m-%d %H:%M:%S" '
                        f'(for "before" and "after" only)') from err
        try:
            if dt < datetime(1900, 1, 1):
                logger.warning('Timestamp is before 1900-JAN-01, ignoring')
                return False
        except (ValueError, OverflowError):
            logger.warning('Could not parse timestamp, ignoring')
            return False
        return True

    ek = {'on', 'before', 'after'}
    if sum(k in ek for k in ts_dict) > 0:
        if 'on' in ts_dict and \
                sum(k in ek for k in ts_dict) > 1 and \
                _is_valid_ts('on', ts_dict['on']):
            logger.warning('Ignoring any other keys than "on"')
            ts = {'on': ts_dict['on']}
        else:
            ts = {k: v for k, v in ts_dict.items() if k in ek and
                  _is_valid_ts(k, v)}
    else:
        raise ValueError(f'None of the allowed keys '
                         f'{", ".join(list(ek))} were provided')
    return ts


def _jsonify_query_data(readers=None, versions=None, document_ids=None,
                        timestamp=None):
    """Check and json.dumps the metadata dictionary

    Parameters
    ----------
    readers : list
        The list of reading systems.
    versions : list
        Versions of reading systems.
    document_ids : list
        Document IDs.
    timestamp : dict("on"|"before"|"after",str)
        Reader output time stamp constraint.

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
        pd['timestamp'] = _check_timestamp_dict(timestamp)
    elif timestamp is not None:
        raise ValueError('Argument "timestamp" must be of type dict')

    return json.dumps(pd)
