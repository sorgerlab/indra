import requests

from indra import get_config
from indra.statements import stmts_from_json


def get_statements(subject=None, object=None, agents=None, stmt_type=None):
    agent_strs = ['agent=%s' % agent_str for agent_str in agents]
    params = {}
    for param_key, param_val in [('subject', subject), ('object', object),
                                 ('type', stmt_type)]:
        if param_val is not None:
            params[param_key] = param_val
    resp = submit_request(*agent_strs, **params)
    stmts_json = resp.json()
    return stmts_from_json(stmts_json)


def submit_request(*args, **kwargs):
    query_str = '?' + '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()]
                               + list(args))
    url = get_config('INDRA_DB_REST_URL', failure_ok=False)
    api_key = get_config('INDRADB_REST_API_KEY', failure_ok=False)
    resp = requests.get(url + '/statements/' + query_str,
                        headers={'x-api-key': api_key})
    return resp