import logging
import requests
from indra.config import get_config


logger = logging.getLogger(__name__)


dart_uname = get_config('DART_WM_USERNAME')
dart_pwd = get_config('DART_WM_PASSWORD')


dart_url = 'https://indra-ingest-pipeline-rest-1.prod.dart.worldmodelers.com' \
           '/dart/api/v1/readers/query'


def query_dart_notifications(readers=None, versions=None, document_ids=None,
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
    """
    if all(v is None for v in [readers, versions, document_ids, timestamp]):
        return {}
    pd = {}
    if readers:
        pd['readers'] = readers
    if versions:
        pd['versions'] = versions
    if document_ids:
        pd['document_ids'] = document_ids
    if isinstance(timestamp, dict):
        pass  # Check
    res = requests.post(
        dart_url,
        data={'metadata':
                None
              },
        auth=(dart_uname, dart_pwd)
    )
    if res.status_code != 200:
        logger.warning(f'Dart Notifications Endpoint returned with status'
                       f' {res.status_code}: {res.text}')

        return {}
    return res.json()
