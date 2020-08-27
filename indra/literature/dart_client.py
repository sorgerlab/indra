"""A client for accessing reader output from the DART system."""

import json
import logging
import requests
import itertools
from datetime import datetime
from collections import defaultdict
from indra.config import get_config


logger = logging.getLogger(__name__)


dart_uname = get_config('DART_WM_USERNAME')
dart_pwd = get_config('DART_WM_PASSWORD')
# The URL is configurable since it is subject to change per use case
dart_base_url = get_config('DART_WM_URL')
if dart_base_url is None:
    dart_base_url = ('https://wm-ingest-pipeline-rest-1.prod.dart'
                     '.worldmodelers.com/dart/api/v1/readers')
meta_endpoint = dart_base_url + '/query'
downl_endpoint = dart_base_url + '/download/'


def get_content_by_storage_key(storage_key):
    """Return content from DART based on its storage key.

    Parameters
    ----------
    storage_key : str
        A DART storage key.

    Returns
    -------
    dict
        The content corresponding to the storage key.
    """
    res = requests.get(url=downl_endpoint + storage_key,
                       auth=(dart_uname, dart_pwd))
    res.raise_for_status()
    return res.text


def get_reader_outputs(readers=None, versions=None, document_ids=None,
                       timestamp=None):
    """Return reader outputs by querying the DART API.

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
    dict(str, dict)
        A two-level dict of reader output keyed by reader and then
        document id.
    """
    records = get_reader_output_records(readers=readers, versions=versions,
                                        document_ids=document_ids,
                                        timestamp=timestamp)
    logger.info('Got %d document storage keys. Fetching output...' %
                len(records))
    return download_records(records)


def download_records(records):
    """Return reader outputs corresponding to a list of records.

    Parameters
    ----------
    records : list of dict
        A list of records returned from the reader output query.

    Returns
    -------
    dict(str, dict)
        A two-level dict of reader output keyed by reader and then
        document id.
    """
    # Loop document keys and get documents
    reader_outputs = defaultdict(dict)
    for record in records:
        reader = record['identity']
        doc_id = record['document_id']
        storage_key = record['storage_key']
        try:
            reader_outputs[reader][doc_id] = \
                get_content_by_storage_key(storage_key)
        except Exception as e:
            logger.warning('Error downloading %s' % storage_key)
    return dict(reader_outputs)


def prioritize_records(records, priorities=None):
    """Return unique records per reader and document prioritizing by version.

    Parameters
    ----------
    records : list of dict
        A list of records returned from the reader output query.
    priorities : dict of list
        A dict keyed by reader names (e.g., cwms, eidos) with values
        representing reader versions in decreasing order of priority.

    Returns
    -------
    records : list of dict
        A list of records that are unique per reader and document, picked by
        version priority when multiple records exist for the same reader
        and document.
    """
    priorities = {} if not priorities else priorities
    prioritized_records = []
    key = lambda x: (x['identity'], x['document_id'])
    for (reader, doc_id), group in itertools.groupby(sorted(records, key=key),
                                                     key=key):
        group_records = list(group)
        if len(group_records) == 1:
            prioritized_records.append(group_records[0])
        else:
            reader_prio = priorities.get(reader)
            if reader_prio:
                first_rec = sorted(
                    group_records,
                    key=lambda x: reader_prio.index(x['version']))[0]
                prioritized_records.append(first_rec)
            else:
                logger.warning('Could not prioritize between records: %s' %
                               str(group_records))
                prioritized_records.append(group_records[0])
    return prioritized_records


def get_reader_output_records(readers=None, versions=None, document_ids=None,
                              timestamp=None):
    """Return reader output metadata records by querying the DART API

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
    full_query_data = {'metadata': query_data}
    res = requests.post(meta_endpoint, data=full_query_data,
                        auth=(dart_uname, dart_pwd))
    res.raise_for_status()
    rj = res.json()

    # This handles both empty list and dict
    if not rj or 'records' not in rj:
        return []
    return rj['records']


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
