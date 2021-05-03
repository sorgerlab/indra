"""
Retrieving the results of large queries from the INDRA Database REST API
generally involves multiple individual calls. The Processor classes
defined here manage the retrieval process for results of two types, Statements
and Statement hashes. Instances of these Processors are returned by the query
functions in :py:mod:`indra.sources.indra_db_rest.api`.
"""

import logging
from copy import deepcopy

from threading import Thread
from datetime import datetime

from requests import Timeout

from indra.statements import stmts_from_json
from indra.util.statement_presentation import get_available_source_counts, \
    get_available_ev_counts

from .query import Query
from .util import RecordableLogger
from .util import logger as util_logger
from .exceptions import IndraDBRestResponseError

logger = logging.getLogger('indra_db_rest.query_processor')
request_logger = RecordableLogger('indra_db_rest.request_logs')


class IndraDBQueryProcessor:
    """The parent of all db query processors.

    Parameters
    ----------
    query : :py:class:`Query`
        The query to be evaluated in return for statements.
    limit : int or None
        Select the maximum number of statements to return. When set less than
        500 the effect is much the same as setting persist to false, and will
        guarantee a faster response. Default is None.
    sort_by : str or None
        Options are currently 'ev_count' or 'belief'. Results will return in
        order of the given parameter. If None, results will be turned in an
        arbitrary order.
    persist : bool
        Default is True. When False, if a query comes back limited (not all
        results returned), just give up and pass along what was returned.
        Otherwise, make further queries to get the rest of the data (which may
        take some time).
    timeout : positive int or None
        If an int, return after `timeout` seconds, even if query is not done.
        Default is None.
    strict_stop : bool
        If True, the query will only be given timeout to complete before being
        abandoned entirely. Otherwise the timeout will simply wait for the
        thread to join for `timeout` seconds before returning, allowing other
        work to continue while the query runs in the background. The default is
        False. NOTE: in practice, due to overhead, the precision of the timeout
        is only around +/-0.1 seconds.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can
        also help gracefully handle an unreliable connection, if you're
        willing to wait. Default is 3
    api_key : str or None
        Override or use in place of the API key given in the INDRA config file.
    """
    result_type = NotImplemented

    def __init__(self, query: Query, limit=None, sort_by='ev_count',
                 timeout=None, strict_stop=False, persist=True, tries=3,
                 api_key=None):
        self.query = query
        self.limit = limit
        self.sort_by = sort_by
        self.tries = tries
        self.__strict_stop = strict_stop
        self.__timeout = timeout
        self.__offset = 0
        self.__quota = limit
        self.__api_key = api_key
        self.__canceled = False
        self.__start_time = None
        self.__th = None
        self.requests_completed = 0

        self._evidence_counts = {}
        self._belief_scores = {}
        self._source_counts = {}

        if limit != 0:
            self._run(persist=persist)

    # Metadata Retrieval methods.

    def get_ev_counts(self):
        """Get a dictionary of evidence counts."""
        return self._evidence_counts.copy()

    def get_belief_scores(self):
        """Get a dictionary of belief scores."""
        return self._belief_scores.copy()

    def get_source_counts(self):
        """Get the source counts as a dict per statement hash."""
        return deepcopy(self._source_counts)

    # Process control methods

    def cancel(self):
        """Cancel the job, stopping the thread running in the background."""
        self.__canceled = True

    def is_working(self):
        """Check if the thread is running."""
        if not self.__th:
            return False
        return self.__th.is_alive()

    def wait_until_done(self, timeout=None):
        """Wait for the background load to complete."""
        if not self.__th:
            raise IndraDBRestResponseError("There is no thread waiting to "
                                           "complete.")
        start = datetime.now()
        self.__th.join(timeout)
        dt = datetime.now() - start
        if self.__th.is_alive():
            logger.warning("Timed out after %0.3f seconds waiting for "
                           "statement load to complete." % dt.total_seconds())
            ret = False
        else:
            logger.info("Waited %0.3f seconds for statements to finish "
                        "loading." % dt.total_seconds())
            ret = True
        return ret

    @staticmethod
    def print_quiet_logs():
        """Print the logs that were suppressed during the query."""
        print(request_logger.get_quiet_logs())

    # Helper methods

    def _get_next_offset(self):
        """Get the offset of the next web request that will be made."""
        return self.__offset

    def _get_next_limit(self):
        """Get the limit of the next web request that will be made."""
        return self.__quota

    def _mark_start(self):
        self.__start_time = datetime.now()

    def _time_since_start(self):
        dt = datetime.now() - self.__start_time
        return dt.total_seconds()

    def _strict_time_is_up(self):
        if self.__start_time is not None and self.__strict_stop:
            if self._time_since_start() > self.__timeout:
                return True
        return False

    def _done(self):
        return (self.__canceled
                or self.__offset is None
                or self.__offset > 0 and self.__quota == 0
                or self._strict_time_is_up())

    def _set_special_params(self, **params):
        self.__special_params = params

    def _run_query(self):
        # If we are in strict stop mode, we want to be sure we give up after
        # the given overall timeout, so we need to account for time spend on
        # other queries.
        if self.__strict_stop:
            query_timeout = self.__timeout - self._time_since_start()
            if query_timeout <= 0:
                return
        else:
            query_timeout = None

        # Run the query.
        try:
            r = self.requests_completed
            nth = f"{r}{['st', 'nd', 'rd'][r-1] if 0 < r < 4 else 'th'}"
            request_logger.info(f"Running {nth} request for {self.result_type}")
            request_logger.info(f"  LIMIT: {self.__quota}")
            request_logger.info(f"  OFFSET: {self.__offset}")
            if query_timeout:
                request_logger.info(f"  TIMEOUT: {query_timeout}")

            result = self.query.get(self.result_type, offset=self.__offset,
                                    limit=self.__quota, sort_by=self.sort_by,
                                    timeout=query_timeout, n_tries=self.tries,
                                    api_key=self.__api_key,
                                    **self.__special_params)
        except Timeout:
            # Make sure this is the timeout we think it is.
            if not self.__strict_stop or not self._strict_time_is_up():
                raise
            return

        # Update results
        self._evidence_counts.update(result.evidence_counts)
        self._belief_scores.update(result.belief_scores)
        self._handle_new_result(result, self._source_counts)

        # Update the quota
        if self.__quota is not None:
            self.__quota -= len(result.results)

        # Increment the page
        self.__offset = result.next_offset

        # Increment the number of queries run.
        self.requests_completed += 1

        return

    def _run_queries(self, persist):
        """Use paging to get all statements requested."""
        self._mark_start()
        self._run_query()

        # Check if we want to keep going.
        if not persist:
            self._compile_results()
            return

        # Get the rest of the content.
        while not self._done():
            self._run_query()

        # Create the actual statements.
        self._compile_results()

        # This is end of the loop, one way or another. Restore logging if it
        # was redirected.
        request_logger.unquiet()
        util_logger.unquiet()
        return

    def _run(self, persist=True):
        util_logger.quiet()
        query_english = self.query.get_query_english()
        logger.info(f"Retrieving {self.result_type} that {query_english}.")

        # Handle the content if we were limited.
        self.__th = Thread(target=self._run_queries,
                           args=[persist])
        self.__th.start()

        if self.__timeout is None:
            logger.debug("Waiting for thread to complete...")
            self.__th.join()
        else:
            if self.__timeout:  # is not 0
                logger.debug("Waiting at most %d seconds for thread to"
                             "complete..." % self.__timeout)
                self.__th.join(self.__timeout)
            if not self._done():
                request_logger.quiet()
                logger.info("Leaving request to background thread. Logs "
                            "may be viewed using the `print_quiet_logs()` "
                            "method.")
        return

    # Child defined methods

    def _compile_results(self):
        raise NotImplementedError()

    def _handle_new_result(self, result, source_counts):
        raise NotImplementedError()


class DBQueryStatementProcessor(IndraDBQueryProcessor):
    """A Processor to get Statements from the server.

    For information on thread control and other methods, see the docs for
    :py:class:`IndraDBQueryProcessor`.

    Parameters
    ----------
    query : :py:class:`Query`
        The query to be evaluated in return for statements.
    limit : int or None
        Select the maximum number of statements to return. When set less than
        500 the effect is much the same as setting persist to false, and will
        guarantee a faster response. Default is None.
    ev_limit : int or None
        Limit the amount of evidence returned per Statement. Default is 100.
    filter_ev : bool
        Indicate whether evidence should have the same filters applied as
        the statements themselves, where appropriate (e.g. in the case of a
        filter by paper).
    sort_by : str or None
        Options are currently 'ev_count' or 'belief'. Results will return in
        order of the given parameter. If None, results will be turned in an
        arbitrary order.
    persist : bool
        Default is True. When False, if a query comes back limited (not all
        results returned), just give up and pass along what was returned.
        Otherwise, make further queries to get the rest of the data (which may
        take some time).
    timeout : positive int or None
        If an int, return after `timeout` seconds, even if query is not done.
        Default is None.
    strict_stop : bool
        If True, the query will only be given timeout to complete before being
        abandoned entirely. Otherwise the timeout will simply wait for the
        thread to join for `timeout` seconds before returning, allowing other
        work to continue while the query runs in the background. The default is
        False.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can
        also help gracefully handle an unreliable connection, if you're
        willing to wait. Default is 3.
    api_key : str or None
        Override or use in place of the API key given in the INDRA config file.

    """
    result_type = 'statements'

    def __init__(self, query: Query, limit=None, sort_by='ev_count',
                 ev_limit=10, filter_ev=True, timeout=None, strict_stop=False,
                 persist=True, use_obtained_counts=False, tries=3,
                 api_key=None):

        self.statements = []
        self.statements_sample = None

        self.__statement_jsons = {}
        self.__started = False

        self.use_obtained_counts = use_obtained_counts
        self._set_special_params(ev_limit=ev_limit, filter_ev=filter_ev)
        super(DBQueryStatementProcessor, self).\
            __init__(query, limit=limit, sort_by=sort_by, timeout=timeout,
                     strict_stop=strict_stop, persist=persist, tries=tries,
                     api_key=api_key)

    # Metadata Retrieval methods.

    def get_ev_count_by_hash(self, stmt_hash):
        """Get the total evidence count for a statement hash."""
        return self._evidence_counts.get(stmt_hash, 0)

    def get_ev_count(self, stmt):
        """Get the total evidence count for a statement."""
        return self.get_ev_count_by_hash(stmt.get_hash(shallow=True))

    def get_belief_score_by_hash(self, stmt_hash):
        """Get the belief score for a statement hash."""
        return self._belief_scores.get(stmt_hash, 0)

    def get_belief_score_by_stmt(self, stmt):
        """Get the belief score for a statement."""
        return self.get_belief_score_by_hash(stmt.get_hash(shallow=True))

    def get_hash_statements_dict(self):
        """Return a dict of Statements keyed by hashes."""
        res = {stmt_hash: stmts_from_json([stmt])[0]
               for stmt_hash, stmt in self.__statement_jsons.items()}
        return res

    def get_source_count_by_hash(self, stmt_hash):
        """Get the source counts for a given statement."""
        return self._source_counts.get(stmt_hash, {})

    def get_source_count(self, stmt):
        """Get the source counts for a given statement."""
        return self.get_source_count_by_hash(stmt.get_hash(shallow=True))

    # Result merging methods

    def merge_results(self, other_processor):
        """Merge the results of this processor with those of another."""
        if not isinstance(other_processor, self.__class__):
            raise ValueError(f"Can only extend with another "
                             f"{self.__class__.__name__} instance.")

        # Where there is overlap, there _should_ be agreement.
        self._evidence_counts.update(other_processor._evidence_counts)
        self._source_counts.update(other_processor._source_counts)
        self._belief_scores.update(other_processor._belief_scores)

        # Merge the statement JSONs.
        for k, sj in other_processor.__statement_jsons.items():
            if k not in self.__statement_jsons:
                self.__statement_jsons[k] = sj  # This should be most of them
            else:
                # This should only happen rarely.
                for evj in sj['evidence']:
                    self.__statement_jsons[k]['evidence'].append(evj)

        # Recompile the statements
        self._compile_results()
        return

    # Helper methods

    def _handle_new_result(self, result, source_counts):
        """Merge these statement jsons with new jsons."""
        # Merge counts.
        source_counts.update(result.source_counts)

        # Merge JSONs
        for k, sj in result.results.items():
            if k not in self.__statement_jsons:
                self.__statement_jsons[k] = sj  # This should be most of them
            else:
                # This should only happen rarely.
                for evj in sj['evidence']:
                    self.__statement_jsons[k]['evidence'].append(evj)

        # Add to the sample.
        if not self.__started:
            self.statements_sample = stmts_from_json(result.results.values())
            self.__started = True
        return

    def _compile_results(self):
        """Generate statements from the jsons."""
        self.statements = stmts_from_json(self.__statement_jsons.values())
        if self.use_obtained_counts:
            self.__source_counts = get_available_source_counts(self.statements)
            self.__evidence_counts = get_available_ev_counts(self.statements)


class DBQueryHashProcessor(IndraDBQueryProcessor):
    """A processor to get hashes from the server.

    Parameters
    ----------
    query : :py:class:`Query`
        The query to be evaluated in return for statements.
    limit : int or None
        Select the maximum number of statements to return. When set less than
        500 the effect is much the same as setting persist to false, and will
        guarantee a faster response. Default is None.
    sort_by : str or None
        Options are currently 'ev_count' or 'belief'. Results will return in
        order of the given parameter. If None, results will be turned in an
        arbitrary order.
    persist : bool
        Default is True. When False, if a query comes back limited (not all
        results returned), just give up and pass along what was returned.
        Otherwise, make further queries to get the rest of the data (which may
        take some time).
    timeout : positive int or None
        If an int, return after `timeout` seconds, even if query is not done.
        Default is None.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can
        also help gracefully handle an unreliable connection, if you're
        willing to wait. Default is 3.
    """
    result_type = 'hashes'

    def __init__(self, *args, **kwargs):
        self.hashes = []
        super(DBQueryHashProcessor, self).__init__(*args, **kwargs)

    def _handle_new_result(self, result, source_counts):
        source_counts.update(result.source_counts)
        self.hashes.extend(result.results)

    def _compile_results(self):
        pass
