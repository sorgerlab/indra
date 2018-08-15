from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_statements', 'get_statements_for_paper',
           'get_statements_by_hash', 'IndraDBRestError']

import json
import logging
import requests
from os.path import join
from threading import Thread

from indra.db.util import _clockit

from indra import get_config
from indra.statements import stmts_from_json, get_statement_by_name, \
    get_all_descendants


logger = logging.getLogger('db_rest_client')


class IndraDBRestError(Exception):
    def __init__(self, resp):
        self.status_code = resp.status_code
        if hasattr(resp, 'text'):
            self.reason = resp.text
        else:
            self.reason = resp.reason
        reason_lines = self.reason.splitlines()
        if len(reason_lines) > 1:
            ws = '\n'
            for i, line in enumerate(reason_lines):
                reason_lines[i] = '  ' + line
        else:
            ws = ' '
        fmtd_reason = '\n'.join(reason_lines)
        self.resp = resp
        Exception.__init__(self, ('Got bad return code %d:%s%s'
                                  % (self.status_code, ws, fmtd_reason)))
        return


class IndraDBRestResponse(object):
    """The packaging for query responses."""
    def __init__(self, statement_jsons):
        self.statements = []
        self.statements_sample = None
        self.statement_jsons = statement_jsons
        self.done = False
        return

    def extend_statements(self, other_response):
        if not isinstance(other_response, self.__class__):
            raise ValueError("Can only extend with another %s instance."
                             % self.__class__.__name__)

        self.statements.extend(other_response.statements)
        if other_response.statements_sample is not None:
            if self.statements_sample is None:
                self.statements_sample = other_response.statements_sample
            else:
                self.statements_sample.extend(other_response.statements_sample)
        return

    def merge_json(self, stmt_json):
        for k, sj in stmt_json.items():
            if k not in self.statement_jsons:
                self.statement_jsons[k] = sj  # This should be most of them
            else:
                # This should only happen rarely.
                for evj in sj['evidence']:
                    self.statement_jsons[k]['evidence'].append(evj)
        return

    def compile_statements(self):
        self.statements = stmts_from_json(self.statement_jsons.values())


def _query_and_extract(agent_strs, params):
    resp = _submit_query_request('statements', *agent_strs, **params)
    resp_dict = resp.json()
    stmts_json = resp_dict['statements']
    total_ev = resp_dict['total_evidence']
    limit_ev = resp_dict['limit']

    # NOTE: this is technically not a direct conclusion, and could be wrong,
    # resulting in an unnecessary extra query, but that should almost never
    # happen.
    limited = (total_ev >= limit_ev)
    return stmts_json, total_ev, limited


def _get_statements_persistently(agent_strs, params, offset, total_ev, ret):
    """Use paging to get all statements."""
    limited = True

    # Get the rest of the content.
    while limited:
        # Get the next page.
        offset = offset + total_ev
        params['offset'] = offset
        new_stmts_json, new_total_ev, limited = \
            _query_and_extract(agent_strs, params)

        # Merge in the new results
        ret.merge_json(new_stmts_json)

    # Create the actual statements.
    ret.compile_statements()
    ret.done = True
    return


def _make_stmts_query(agent_strs, params, persist=True, block=True):
    """Slightly lower level function to get statements from the REST API."""
    # Perform the first (and last?) query
    stmts_json, total_ev, limited = _query_and_extract(agent_strs, params)

    # Initialize the return dict.
    ret = IndraDBRestResponse(stmts_json)

    # Handle the content if we were limited.
    if limited:
        logger.info("Some results could not be returned directly.")
        if persist:
            offset = params.get('offset', 0)
            args = [agent_strs, params, offset, total_ev, ret]
            if block:
                logger.info("You chose to persist, so I will paginate through "
                            "the rest until I have everything!")
                _get_statements_persistently(*args)
                assert ret.done
            else:
                logger.info("You chose to persist without blocking. Pagination "
                            "is being performed in a thread.")
                ret.statements_sample = stmts_from_json(stmts_json.values())
                th = Thread(target=_get_statements_persistently, args=args)
                th.start()
        else:
            logger.warning("You did not choose persist=True, therefore this is "
                           "all you get.")
            ret.compile_statements()
            ret.statements_sample = ret.statements[:]
            ret.done = True
    else:
        ret.compile_statements()
        ret.done = True
    return ret


@_clockit
def get_statements(subject=None, object=None, agents=None, stmt_type=None,
                   use_exact_type=False, offset=None, persist=True, block=True,
                   simple_response=True):
    """Get statements from INDRA's database using the web api.

    Parameters
    ----------
    subject/object : str
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
        sensitive. If the statement class given has sub-classes
        (e.g. RegulateAmount has IncreaseAmount and DecreaseAmount), then both
        the class itself, and its subclasses, will be queried, by default. If
        you do not want this behavior, set use_exact_type=True.
    use_exact_type : bool
        If stmt_type is given, and you only want to search for that specific
        statement type, set this to True. Default is False.
    on_limit : str
        There are four options for handling a query that is too large:
        `sample` - (default) take a sample of statements from the result,
        `truncate` - simply return the first 10,000 statements of the result,
        `error` - raise an error if the query is too large, or
        `persist` - perform as many queries as needed to get all the statements.
        Note that this last option generally takes much much longer to execute.

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statement instances. Note that if a supporting or
        supported Statement was not included in your query, it will simply be
        instantiated as an `Unresolved` statement, with `uuid` of the statement.
    """
    # Make sure we got at least SOME agents (the remote API will error if we
    # we proceed with no arguments.
    if subject is None and object is None and agents is None:
        raise ValueError("At least one agent must be specified, or else "
                         "the scope will be too large.")

    # Formulate inputs for the agents..
    agent_strs = [] if agents is None else ['agent=%s' % ag for ag in agents]
    key_val_list = [('subject', subject), ('object', object)]
    params = {param_key: param_val for param_key, param_val in key_val_list
              if param_val is not None}
    if offset is not None:
        params['offset'] = offset

    # Handle the type(s).
    if stmt_type is not None:
        if use_exact_type:
            params['type'] = stmt_type
            resp = _make_stmts_query(agent_strs, params, persist, block)
        else:
            stmt_class = get_statement_by_name(stmt_type)
            descendant_classes = get_all_descendants(stmt_class)
            stmt_types = [cls.__name__ for cls in descendant_classes] \
                + [stmt_type]
            resp = None
            for stmt_type in stmt_types:
                params['type'] = stmt_type
                new_resp = _make_stmts_query(agent_strs, params, persist)
                logger.info("Found %d %s statements."
                            % (len(resp['statements']), stmt_type))
                if resp is None:
                    resp = new_resp
                else:
                    resp.extend_statements(new_resp)
    else:
        resp = _make_stmts_query(agent_strs, params, persist, block)

    if simple_response:
        ret = resp.statements
    else:
        ret = resp
    return ret


@_clockit
def get_statements_by_hash(hash_list, with_support=True):
    """Get fully formed statements from a list of hashes."""
    resp = _submit_request('post', 'statements/from_hashes',
                           hashes=hash_list)
    return stmts_from_json(resp.json()['statements'].values())


@_clockit
def get_statements_for_paper(id_val, id_type='pmid'):
    """Get the set of raw Statements extracted from a paper given by the id.

    Note that in the future this will draw upon preassembled Statements.

    Parameters
    ----------
    id_val : str or int
        The value of the id of the paper of interest.
    id_type : str
        This may be one of 'pmid', 'pmcid', 'doi', 'pii', 'manuscript id', or
        'trid', which is the primary key id of the text references in the
        database. The default is 'pmid'.

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statement instances.
    """
    resp = _submit_query_request('papers', id=id_val, type=id_type)
    stmts_json = resp.json()['statements']
    return stmts_from_json(stmts_json)


def _submit_query_request(end_point, *args, **kwargs):
    """Low level function to format the query string."""
    # This isn't handled by requests because of the multiple identical agent
    # keys, e.g. {'agent': 'MEK', 'agent': 'ERK'} which is not supported in
    # python, but is allowed and necessary in these query strings.
    query_str = '?' + '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()]
                               + list(args))
    return _submit_request('get', end_point, query_str)


@_clockit
def _submit_request(meth, end_point, query_str='', **data):
    """Even lower level function to make the request."""
    url = get_config('INDRA_DB_REST_URL', failure_ok=False)
    api_key = get_config('INDRA_DB_REST_API_KEY', failure_ok=False)
    url_path = join(url, end_point)
    if query_str:
        url_path += '/' + query_str
    headers={'x-api-key': api_key}
    if data:
        # This is an assumption which applies to our use cases for now, but may
        # not generalize.
        headers['content-type'] = 'application/json'
        json_data = json.dumps(data)
    else:
        json_data = None
    print('url and query string:', url_path)
    print('headers:', headers)
    print('data:', data)
    method_func = getattr(requests, meth.lower())
    resp = method_func(url_path, headers=headers, data=json_data)
    if resp.status_code == 200:
        return resp
    else:
        raise IndraDBRestError(resp)
