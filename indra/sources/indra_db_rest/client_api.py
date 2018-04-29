from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_statements', 'IndraDBRestError']

import requests

from indra import get_config
from indra.statements import stmts_from_json


class IndraDBRestError(Exception):
    pass


def get_statements(subject=None, object=None, agents=None, stmt_type=None):
    """Get statements from INDRA's database using the web api.

    This tool is intended for those that wish to use the cumulative database
    of pre-assembled INDRA Statements compiled from all our input databases and
    from reading all of the available medical literature, but who do not have
    direct access to the entire database.

    Such access will require you have both a URL (`INDRA_DB_REST_URL`) and
    possibly an API key (`INDRA_DB_REST_API_KEY`), both of which may be placed in
    your config file or as environment variables.

    If you do not have these, but would like to access the database rest api,
    you may contact the developers to request a URL and key.

    Parameters
    ----------
    subject, object : str
        Optionally specify the subject and/or object of the statements in
        you wish to get from the database. By default, the namespace is assumed
        to be HGNC gene names, however you may specify another namespace by
        including `@<namespace>` at the end of the name string. For example, if
        you want to specify an agent by chebi, you could use `CHEBI:6801@CHEBI`,
        or if you wanted to use the HGNC id, you could use `6871@HGNC`.
    agents : list[str]
        A list of agents, specified in the same manner as subject and object,
        but without specifying their grammatical position.
    stmt_type : str
        Specify the types of interactions you are interested in, as indicated
        by the sub-classes of INDRA's Statements. This argument is *not* case
        sensitive.

    Returns
    -------
    stmts : list[:pyclass:`indra.statements.Statement`]
        A list of INDRA Statement instances. Note that if a supporting or
        supported Statement was not included in your query, it will simply be
        instantiated as an `Unresolved` statement, with `uuid` of the statement.
    """
    if subject is None and object is None and agents is None:
        raise IndraDBRestError("At least one agent must be specified, or else "
                               "the scope will be too large.")
    if agents is not None:
        agent_strs = ['agent=%s' % agent_str for agent_str in agents]
    else:
        agent_strs = []
    params = {}
    for param_key, param_val in [('subject', subject), ('object', object),
                                 ('type', stmt_type)]:
        if param_val is not None:
            params[param_key] = param_val
    resp = _submit_request(*agent_strs, **params)
    stmts_json = resp.json()
    return stmts_from_json(stmts_json)


def _submit_request(*args, **kwargs):
    """Low level function to make the request to the rest API."""
    query_str = '?' + '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()]
                               + list(args))
    url = get_config('INDRA_DB_REST_URL', failure_ok=False)
    api_key = get_config('INDRA_DB_REST_API_KEY', failure_ok=False)
    resp = requests.get(url + '/statements/' + query_str,
                        headers={'x-api-key': api_key})
    if resp.status_code == 200:
        return resp
    else:
        raise IndraDBRestError("Got bad return code %d: %s"
                               % (resp.status_code, resp.reason))
