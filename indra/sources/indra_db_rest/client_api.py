from __future__ import absolute_import, unicode_literals
from builtins import dict, str

__all__ = ['get_statements', 'get_statements_for_paper',
           'get_statements_by_hash', 'IndraDBRestAPIError',
           'IndraDBRestClientError', 'IndraDBRestResponseError']

import json
import logging
import requests
from threading import Thread
from datetime import datetime
from collections import OrderedDict, defaultdict
from urllib.parse import urljoin

from indra.util import clockit

from indra import get_config
from indra.statements import stmts_from_json, get_statement_by_name, \
    get_all_descendants


logger = logging.getLogger('db_rest_client')


class IndraDBRestClientError(Exception):
    pass


class IndraDBRestResponseError(IndraDBRestClientError):
    pass


class IndraDBRestAPIError(IndraDBRestClientError):
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
    def __init__(self, max_stmts=None):
        self.statements = []
        self.statements_sample = None
        self.__statement_jsons = {}
        self.__done_dict = defaultdict(lambda: False)
        self.__evidence_counts = {}
        self.__started = False
        self.__page_dict = defaultdict(lambda: 0)
        self.__th = None
        self.__quota = max_stmts
        return

    def is_working(self):
        """Check if the thread is running."""
        if not self.__th:
            return False
        return self.__th.is_alive()

    def get_ev_count(self, stmt):
        """Get the total evidence count for a statement."""
        return self.__evidence_counts.get(str(stmt.get_hash(shallow=True)))

    def get_hash_statements_dict(self):
        """Return a dict of Statements keyed by hashes."""
        res = {stmt_hash: stmts_from_json([stmt])[0]
               for stmt_hash, stmt in self.__statement_jsons.items()}
        return res

    def extend_statements(self, other_response):
        """Extend this object with new statements."""
        if not isinstance(other_response, self.__class__):
            raise ValueError("Can only extend with another %s instance."
                             % self.__class__.__name__)
        self.statements.extend(other_response.statements)
        if other_response.statements_sample is not None:
            if self.statements_sample is None:
                self.statements_sample = other_response.statements_sample
            else:
                self.statements_sample.extend(other_response.statements_sample)

        self.merge_json(other_response.__statement_jsons,
                        other_response.__evidence_counts)
        return

    def merge_json(self, stmt_json, ev_counts):
        """Merge these statement jsons with new jsons."""
        # Where there is overlap, there _should_ be agreement.
        self.__evidence_counts.update(ev_counts)

        for k, sj in stmt_json.items():
            if k not in self.__statement_jsons:
                self.__statement_jsons[k] = sj  # This should be most of them
            else:
                # This should only happen rarely.
                for evj in sj['evidence']:
                    self.__statement_jsons[k]['evidence'].append(evj)

        if not self.__started:
            self.statements_sample = stmts_from_json(
                self.__statement_jsons.values())
            self.__started = True
        return

    def compile_statements(self):
        """Generate statements from the jsons."""
        self.statements = stmts_from_json(self.__statement_jsons.values())

    def wait_until_done(self, timeout=None):
        """Wait for the background load to complete."""
        start = datetime.now()
        if not self.__th:
            raise IndraDBRestResponseError("There is no thread waiting to "
                                           "complete.")
        self.__th.join(timeout)
        now = datetime.now()
        dt = now - start
        if self.__th.is_alive():
            logger.warning("Timed out after %0.3f seconds waiting for "
                           "statement load to complete." % dt.total_seconds())
            ret = False
        else:
            logger.info("Waited %0.3f seconds for statements to finish loading."
                        % dt.total_seconds())
            ret = True
        return ret

    def _all_done(self):
        every_type_done = (len(self.__done_dict) > 0
                           and all(self.__done_dict.values()))
        quota_done = (self.__quota is not None and self.__quota <= 0)
        return every_type_done or quota_done

    def _query_and_extract(self, agent_strs, params, stmt_type=None):
        assert not self._all_done(), "Tried to run query but I'm done!"
        params['offset'] = self.__page_dict[stmt_type]
        params['max_stmts'] = self.__quota
        if stmt_type is not None:
            params['type'] = stmt_type
        resp = _submit_query_request('statements', *agent_strs, **params)
        resp_dict = resp.json(object_pairs_hook=OrderedDict)
        stmts_json = resp_dict['statements']
        ev_totals = resp_dict['evidence_totals']
        page_step = resp_dict['statement_limit']
        num_returned = len(stmts_json)

        # Update the result
        self.merge_json(stmts_json, ev_totals)

        # NOTE: this is technically not a direct conclusion, and could be
        # wrong, resulting in a single unnecessary extra query, but that
        # should almost never happen, and if it does, it isn't the end of
        # the world.
        self.__done_dict[stmt_type] = num_returned < page_step

        # Update the quota
        if self.__quota is not None:
            self.__quota -= num_returned

        # Increment the page
        self.__page_dict[stmt_type] += page_step

        return

    def _query_over_statement_types(self, agent_strs, stmt_types, params):
        if not stmt_types:
            self._query_and_extract(agent_strs, params.copy())
        else:
            for stmt_type in stmt_types:
                if self.__done_dict[stmt_type]:
                    continue
                self._query_and_extract(agent_strs, params.copy(), stmt_type)

                # Check the quota
                if self.__quota is not None and self.__quota <= 0:
                    break
        return

    def _run_queries(self, agent_strs, stmt_types, params, persist):
        """Use paging to get all statements requested."""
        self._query_over_statement_types(agent_strs, stmt_types, params)

        assert len(self.__done_dict) == len(stmt_types) \
            or None in self.__done_dict.keys(), \
            "Done dict was not initiated for all stmt_type's."

        # Check if we want to keep going.
        if not persist:
            self.compile_statements()
            return

        # Get the rest of the content.
        while not self._all_done():
            self._query_over_statement_types(agent_strs, stmt_types, params)

        # Create the actual statements.
        self.compile_statements()
        return

    def make_stmts_queries(self, agent_strs, stmt_types, params, persist=True,
                           block_secs=None):
        """Slightly lower level function to get statements from the REST API."""
        # Handle the content if we were limited.
        logger.info("Some results could not be returned directly.")
        args = [agent_strs, stmt_types, params, persist]
        logger.info("You chose to persist without blocking. Pagination "
                    "is being performed in a thread.")
        self.__th = Thread(target=self._run_queries, args=args)
        self.__th.start()

        if block_secs is None:
            self.__th.join()
        elif block_secs:  # is not 0
            logger.info("Waiting for %d seconds..." % block_secs)
            self.__th.join(block_secs)
        return


@clockit
def get_statements(subject=None, object=None, agents=None, stmt_type=None,
                   use_exact_type=False, persist=True, timeout=None,
                   simple_response=True, ev_limit=10, best_first=True, tries=2,
                   max_stmts=None):
    """Get Statements from the INDRA DB web API matching given agents and type.

    There are two types of response available. You can just get a list of
    INDRA Statements, or you can get an IndraRestResponse object, which allows
    Statements to be loaded in a background thread, providing a sample of the
    best* content available promptly in the sample_statements attribute, and
    populates the statements attribute when the paged load is complete.

    *In the sense of having the most supporting evidence.

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
        you do not want this behavior, set use_exact_type=True. Note that if
        max_stmts is set, it is possible only the exact statement type will
        be returned, as this is the first searched. The processor then cycles
        through the types, getting a page of results for each type and adding it
        to the quota, until the max number of statements is reached.
    use_exact_type : bool
        If stmt_type is given, and you only want to search for that specific
        statement type, set this to True. Default is False.
    persist : bool
        Default is True. When False, if a query comes back limited (not all
        results returned), just give up and pass along what was returned.
        Otherwise, make further queries to get the rest of the data (which may
        take some time).
    timeout : positive int or None
        If an int, block until the work is done and statements are retrieved, or
        until the timeout has expired, in which case the results so far will be
        returned in the response object, and further results will be added in
        a separate thread as they become available. If simple_response is True,
        all statements available will be returned. Otherwise (if None), block
        indefinitely until all statements are retrieved. Default is None.
    simple_response : bool
        If True, a simple list of statements is returned (thus block should also
        be True). If block is False, only the original sample will be returned
        (as though persist was False), until the statements are done loading, in
        which case the rest should appear in the list. This behavior is not
        encouraged. Default is True (for the sake of backwards compatibility).
    ev_limit : int or None
        Limit the amount of evidence returned per Statement. Default is 10.
    best_first : bool
        If True, the preassembled statements will be sorted by the amount of
        evidence they have, and those with the most evidence will be
        prioritized. When using `max_stmts`, this means you will get the "best"
        statements. If False, statements will be queried in arbitrary order.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can also
        help gracefully handle an unreliable connection, if you're willing to
        wait. Default is 2.
    max_stmts : int or None
        Select the maximum number of statements to return. When set less than
        1000 the effect is much the same as setting persist to false, and will
        guarantee a faster response. Default is None.

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
    agent_strs = [] if agents is None else ['agent%d=%s' % (i, ag)
                                            for i, ag in enumerate(agents)]
    key_val_list = [('subject', subject), ('object', object)]
    params = {param_key: param_val for param_key, param_val in key_val_list
              if param_val is not None}
    params['best_first'] = best_first
    params['ev_limit'] = ev_limit
    params['tries'] = tries

    # Handle the type(s).
    stmt_types = [stmt_type] if stmt_type else []
    if stmt_type is not None and not use_exact_type:
        stmt_class = get_statement_by_name(stmt_type)
        descendant_classes = get_all_descendants(stmt_class)
        stmt_types += [cls.__name__ for cls in descendant_classes]

    # Get the response object
    resp = IndraDBRestResponse(max_stmts=max_stmts)
    resp.make_stmts_queries(agent_strs, stmt_types, params, persist, timeout)

    # Format the result appropriately.
    if simple_response:
        ret = resp.statements
    else:
        ret = resp
    return ret


@clockit
def get_statements_by_hash(hash_list, ev_limit=100, best_first=True, tries=2):
    """Get fully formed statements from a list of hashes.

    Parameters
    ----------
    hash_list : list[int or str]
        A list of statement hashes.
    ev_limit : int or None
        Limit the amount of evidence returned per Statement. Default is 100.
    best_first : bool
        If True, the preassembled statements will be sorted by the amount of
        evidence they have, and those with the most evidence will be
        prioritized. When using `max_stmts`, this means you will get the "best"
        statements. If False, statements will be queried in arbitrary order.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can also
        help gracefully handle an unreliable connection, if you're willing to
        wait. Default is 2.
    """
    if not isinstance(hash_list, list):
        raise ValueError("The `hash_list` input is a list, not %s."
                         % type(hash_list))
    if not hash_list:
        return []
    if isinstance(hash_list[0], str):
        hash_list = [int(h) for h in hash_list]
    if not all([isinstance(h, int) for h in hash_list]):
        raise ValueError("Hashes must be ints or strings that can be converted "
                         "into ints.")
    resp = _submit_request('post', 'statements/from_hashes',
                           data={'hashes': hash_list}, ev_limit=ev_limit,
                           best_first=best_first, tries=tries, div='')
    return stmts_from_json(resp.json()['statements'].values())


@clockit
def get_statements_for_paper(id_val, id_type='pmid', ev_limit=10,
                             best_first=True, tries=2, max_stmts=None):
    """Get the set of raw Statements extracted from a paper given by the id.

    Parameters
    ----------
    id_val : str or int
        The value of the id of the paper of interest.
    id_type : str
        This may be one of 'pmid', 'pmcid', 'doi', 'pii', 'manuscript id', or
        'trid', which is the primary key id of the text references in the
        database. The default is 'pmid'.
    ev_limit : int or None
        Limit the amount of evidence returned per Statement. Default is 10.
    best_first : bool
        If True, the preassembled statements will be sorted by the amount of
        evidence they have, and those with the most evidence will be
        prioritized. When using `max_stmts`, this means you will get the "best"
        statements. If False, statements will be queried in arbitrary order.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can also
        help gracefully handle an unreliable connection, if you're willing to
        wait. Default is 2.
    max_stmts : int or None
        Select a maximum number of statements to be returned. Default is None.

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statement instances.
    """
    resp = _submit_query_request('papers', id=id_val, type=id_type,
                                 ev_limit=ev_limit, best_first=best_first,
                                 tries=tries, max_stmts=max_stmts)
    stmts_json = resp.json()['statements']
    return stmts_from_json(stmts_json.values())


def _submit_query_request(end_point, *args, **kwargs):
    """Low level function to format the query string."""
    ev_limit = kwargs.pop('ev_limit', 10)
    best_first = kwargs.pop('best_first', True)
    tries = kwargs.pop('tries', 2)
    # This isn't handled by requests because of the multiple identical agent
    # keys, e.g. {'agent': 'MEK', 'agent': 'ERK'} which is not supported in
    # python, but is allowed and necessary in these query strings.
    query_str = '?' + '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()
                                if v is not None]
                               + list(args))
    return _submit_request('get', end_point, query_str, ev_limit=ev_limit,
                           best_first=best_first, tries=tries)


@clockit
def _submit_request(meth, end_point, query_str='', data=None, ev_limit=50,
                    best_first=True, tries=2, div='/'):
    """Even lower level function to make the request."""
    url = get_config('INDRA_DB_REST_URL', failure_ok=False)
    api_key = get_config('INDRA_DB_REST_API_KEY', failure_ok=True)
    url_path = urljoin(url, end_point)
    if query_str:
        query_str += '&api-key=%s' % api_key
    else:
        query_str = '?api-key=%s' % api_key
    url_path += div + query_str
    headers = {}
    if data:
        # This is an assumption which applies to our use cases for now, but may
        # not generalize.
        headers['content-type'] = 'application/json'
        json_data = json.dumps(data)
    else:
        json_data = None
    logger.info('url and query string: %s',
                url_path.replace(str(api_key), '[api-key]'))
    logger.info('headers: %s', str(headers).replace(str(api_key), '[api-key]'))
    logger.info('data: %s', str(data).replace(str(api_key), '[api-key]'))
    method_func = getattr(requests, meth.lower())
    while tries > 0:
        tries -= 1
        resp = method_func(url_path, headers=headers, data=json_data,
                           params={'ev_limit': ev_limit, 'best_first': best_first})
        if resp.status_code == 200:
            return resp
        elif resp.status_code == 504 and tries > 0:
            logger.warning("Endpoint timed out. Trying again...")
        else:
            raise IndraDBRestAPIError(resp)
