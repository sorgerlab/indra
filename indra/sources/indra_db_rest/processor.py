__all__ = ['IndraDBRestSearchProcessor', 'IndraDBRestHashProcessor']

import logging
from copy import deepcopy
from threading import Thread
from datetime import datetime
from collections import OrderedDict, defaultdict

from indra.statements import stmts_from_json, get_statement_by_name, \
    get_all_descendants

from indra.sources.indra_db_rest.util import submit_query_request, \
    submit_statement_request
from indra.sources.indra_db_rest.exceptions import IndraDBRestResponseError

logger = logging.getLogger(__name__)


class RemoveParam(object):
    pass


class IndraDBRestProcessor(object):
    """The generalized packaging for query responses.

    General Parameters
    ------------------
    timeout : positive int or None
        If an int, block until the work is done and statements are retrieved, or
        until the timeout has expired, in which case the results so far will be
        returned in the response object, and further results will be added in
        a separate thread as they become available. If simple_response is True,
        all statements available will be returned. Otherwise (if None), block
        indefinitely until all statements are retrieved. Default is None.
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

    Attributes
    ----------
    statements : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statements that will be filled once all queries have
        been completed.
    """
    _override_default_api_params = {}

    def __init__(self, *args, **kwargs):
        self.statements = []
        self.__statement_jsons = {}
        self.__evidence_counts = {}
        self.__source_counts = {}

        # Define the basic generic defaults.
        default_api_params = dict(timeout=None, ev_limit=10, best_first=True,
                                  tries=2, max_stmts=None)

        # Update with any overrides.
        default_api_params.update(self._override_default_api_params)

        # Some overrides may be RemoveParam objects, indicating the key should
        # be removed. Filter those out.
        default_api_params = {k: v for k, v in default_api_params.items()
                              if not isinstance(v, RemoveParam)}

        # Update the kwargs to include these default values, if not already
        # specified by the user.
        kwargs.update((k, kwargs.get(k, default_api_params[k]))
                      for k in default_api_params.keys())

        self._run(*args, **kwargs)
        return

    def get_ev_count(self, stmt):
        """Get the total evidence count for a statement."""
        return self.get_ev_count_by_hash(stmt.get_hash(shallow=True))

    def get_ev_count_by_hash(self, stmt_hash):
        """Get the total evidence count for a statement hash."""
        return self.__evidence_counts.get(str(stmt_hash))

    def get_source_counts(self):
        """Get the source counts as a dict per statement hash."""
        return deepcopy(self.__source_counts)

    def get_source_count(self, stmt):
        """Get the source counts for a given statement."""
        return self.get_source_count_by_hash(stmt.get_hash(shallow=True))

    def get_source_count_by_hash(self, stmt_hash):
        """Get the source counts for a given statement."""
        return self.__source_counts.get(stmt_hash)

    def get_ev_counts(self):
        """Get a dictionary of evidence counts."""
        return self.__evidence_counts.copy()

    def get_hash_statements_dict(self):
        """Return a dict of Statements keyed by hashes."""
        res = {stmt_hash: stmts_from_json([stmt])[0]
               for stmt_hash, stmt in self.__statement_jsons.items()}
        return res

    def merge_results(self, other_processor):
        """Merge the results of this processor with those of another."""
        if not isinstance(other_processor, self.__class__):
            raise ValueError("Can only extend with another %s instance."
                             % self.__class__.__name__)
        self.statements.extend(other_processor.statements)
        if other_processor.statements_sample is not None:
            if self.statements_sample is None:
                self.statements_sample = other_processor.statements_sample
            else:
                self.statements_sample.extend(other_processor.statements_sample)

        self._merge_json(other_processor.__statement_jsons,
                         other_processor.__evidence_counts,
                         other_processor.__source_counts)
        return

    def _merge_json(self, stmt_json, ev_counts, source_counts):
        """Merge these statement jsons with new jsons."""
        # Where there is overlap, there _should_ be agreement.
        self.__evidence_counts.update(ev_counts)
        # We turn source counts into an int-keyed dict and update it that way
        self.__source_counts.update({int(k): v
                                     for k, v in source_counts.items()})

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

    def _compile_statements(self):
        """Generate statements from the jsons."""
        self.statements = stmts_from_json(self.__statement_jsons.values())

    def _unload_and_merge_resp(self, resp):
        resp_dict = resp.json(object_pairs_hook=OrderedDict)
        stmts_json = resp_dict['statements']
        ev_totals = resp_dict['evidence_totals']
        source_counts = resp_dict['source_counts']
        limit = resp_dict['statement_limit']
        num_returned = len(stmts_json)

        # Update the result
        self._merge_json(stmts_json, ev_totals, source_counts)

        return limit, num_returned

    def _run(self, *args, **kwargs):
        raise NotImplementedError("_run must be defined in subclass.")


class IndraDBRestHashProcessor(IndraDBRestProcessor):
    """The packaging and processor for hash lookup of statements.

    Parameters
    ----------
    hash_list : list[int or str]
        A list of the matches-key hashes for the statements you want to get.

    Keyword Parameters
    ------------------
    timeout : positive int or None
        If an int, block until the work is done and statements are retrieved, or
        until the timeout has expired, in which case the results so far will be
        returned in the response object, and further results will be added in
        a separate thread as they become available. If simple_response is True,
        all statements available will be returned. Otherwise (if None), block
        indefinitely until all statements are retrieved. Default is None.
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

    Attributes
    ----------
    statements : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statements that will be filled once all queries have
        been completed.
    """
    _default_api_params = {'ev_limit': 100, 'max_stmts': RemoveParam()}

    def _run(self, hash_list, **api_params):
        # Make sure the input is a list (not just a single hash).
        if not isinstance(hash_list, list):
            raise ValueError("The `hash_list` input is a list, not %s."
                             % type(hash_list))

        # If there is nothing in the list, don't waste time with a query.
        if not hash_list:
            return

        # Regularize and check the types of elements in the hash list.
        if isinstance(hash_list[0], str):
            hash_list = [int(h) for h in hash_list]
        if not all([isinstance(h, int) for h in hash_list]):
            raise ValueError("Hashes must be ints or strings that can be "
                             "converted into ints.")

        # Execute the query and load the results.
        resp = submit_statement_request('post', 'from_hashes',
                                        data={'hashes': hash_list},
                                        **api_params)
        self._unload_and_merge_resp(resp)
        self._compile_statements()
        return


class IndraDBRestPaperProcessor(IndraDBRestProcessor):
    """The packaging and processor for hash lookup of statements.

    Parameters
    ----------
    hash_list : list[int or str]
        A list of the matches-key hashes for the statements you want to get.

    Keyword Parameters
    ------------------
    timeout : positive int or None
        If an int, block until the work is done and statements are retrieved, or
        until the timeout has expired, in which case the results so far will be
        returned in the response object, and further results will be added in
        a separate thread as they become available. If simple_response is True,
        all statements available will be returned. Otherwise (if None), block
        indefinitely until all statements are retrieved. Default is None.
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
    max_stmts : int or None
        Select the maximum number of statements to return. When set less than
        1000 the effect is much the same as setting persist to false, and will
        guarantee a faster response. Default is None.

    Attributes
    ----------
    statements : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statements that will be filled once all queries have
        been completed.
    """
    def _run(self, ids, **api_params):
        id_l = [{'id': id_val, 'type': id_type} for id_type, id_val in ids]
        resp = submit_statement_request('post', 'from_papers',
                                        data={'ids': id_l},
                                        **api_params)
        self._unload_and_merge_resp(resp)
        self._compile_statements()
        return


class IndraDBRestSearchProcessor(IndraDBRestProcessor):
    """The packaging for agent and statement type search query responses.

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

    Keyword Parameters
    ------------------
    timeout : positive int or None
        If an int, block until the work is done and statements are retrieved, or
        until the timeout has expired, in which case the results so far will be
        returned in the response object, and further results will be added in
        a separate thread as they become available. If simple_response is True,
        all statements available will be returned. Otherwise (if None), block
        indefinitely until all statements are retrieved. Default is None.
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

    Attributes
    ----------
    statements : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statements that will be filled once all queries have
        been completed.
    statements_sample : list[:py:class:`indra.statements.Statement`]
        A list of the INDRA Statements received from the first query. In
        general these will be the "best" (currently this means they have the
        most evidence) Statements available.
    """
    def is_working(self):
        """Check if the thread is running."""
        if not self.__th:
            return False
        return self.__th.is_alive()

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
            logger.info("Waited %0.3f seconds for statements to finish"
                        "loading." % dt.total_seconds())
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
        resp = submit_query_request('from_agents', *agent_strs, **params)
        page_step, num_returned = self._unload_and_merge_resp(resp)

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
            self._compile_statements()
            return

        # Get the rest of the content.
        while not self._all_done():
            self._query_over_statement_types(agent_strs, stmt_types, params)

        # Create the actual statements.
        self._compile_statements()
        return

    def _run(self, subject=None, object=None, agents=None, stmt_type=None,
             use_exact_type=False, persist=True, **api_params):
        self.statements_sample = None
        self.__started = False
        self.__done_dict = defaultdict(lambda: False)
        self.__page_dict = defaultdict(lambda: 0)
        self.__th = None
        self.__quota = api_params['max_stmts']

        # Make sure we got at least SOME agents (the remote API will error if
        # we proceed with no arguments).
        if subject is None and object is None and not agents:
            raise ValueError("At least one agent must be specified, or else "
                             "the scope will be too large.")

        # Formulate inputs for the agents..
        key_val_list = [('subject', subject), ('object', object)]
        params = {param_key: param_val for param_key, param_val in key_val_list
                  if param_val is not None}
        params.update(api_params)

        agent_strs = [] if agents is None else ['agent%d=%s' % (i, ag)
                                                for i, ag in enumerate(agents)]

        # Handle the type(s).
        stmt_types = [stmt_type] if stmt_type else []
        if stmt_type is not None and not use_exact_type:
            stmt_class = get_statement_by_name(stmt_type)
            descendant_classes = get_all_descendants(stmt_class)
            stmt_types += [cls.__name__ for cls in descendant_classes]

        # Handle the content if we were limited.
        args = [agent_strs, stmt_types, params, persist]
        logger.debug("The remainder of the query will be performed in a "
                     "thread...")
        self.__th = Thread(target=self._run_queries, args=args)
        self.__th.start()

        if api_params['timeout'] is None:
            logger.debug("Waiting for thread to complete...")
            self.__th.join()
        elif api_params['timeout']:  # is not 0
            logger.debug("Waiting at most %d seconds for thread to complete..."
                         % api_params['timeout'])
            self.__th.join(api_params['timeout'])
        return
