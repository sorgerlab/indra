import logging

from threading import Thread
from datetime import datetime

from indra.statements import stmts_from_json
from indra.util.statement_presentation import get_available_source_counts, \
    get_available_ev_counts, standardize_counts

from .query import Query
from .exceptions import IndraDBRestResponseError

logger = logging.getLogger('indra_db_rest.query_processor')


class Processor:
    result_type = NotImplemented

    def __init__(self, query: Query, limit=None, sort_by='default',
                 timeout=None, strict_stop=False, persist=True):
        self.query = query
        self.limit = limit
        self.sort_by = sort_by
        self.__offset = 0
        self.__quota = limit

        self._run(persist=persist, strict_stop=strict_stop, timeout=timeout)

    def _set_special_params(self, **params):
        self.__special_params = params

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

    def _run_query(self, timeout):
        result = self.query.get(self.result_type, offset=self.__offset,
                                limit=self.__quota, sort_by=self.sort_by,
                                timeout=timeout, **self.__special_params)
        self._handle_new_result(result)
        self.__done = result.next_offset is None

        # Update the quota
        if self.__quota is not None:
            self.__quota -= len(result.results)

        # Increment the page
        self.__offset = result.next_offset

        return

    def _run_queries(self, persist, timeout):
        """Use paging to get all statements requested."""
        self._run_query(timeout)

        # Check if we want to keep going.
        if not persist:
            self._compile_results()
            return

        # Get the rest of the content.
        while not self.__done:
            self._run_query(timeout)

        # Create the actual statements.
        self._compile_results()
        return

    def _compile_results(self):
        raise NotImplementedError()

    def _handle_new_result(self, result):
        raise NotImplementedError()

    def _run(self, persist=True, strict_stop=False, timeout=None):
        self.__started = False
        self.__th = None

        # Handle the content if we were limited.
        self.__th = Thread(target=self._run_queries,
                           args=[persist, timeout if strict_stop else None])
        self.__th.start()

        if timeout is None:
            logger.debug("Waiting for thread to complete...")
            self.__th.join()
        elif timeout:  # is not 0
            logger.debug("Waiting at most %d seconds for thread to complete..."
                         % timeout)
            self.__th.join(timeout)
        return


class HashProcessor(Processor):
    """A processor to get hashes from the server."""

    def __init__(self, *args, **kwargs):
        self.hashes = []
        super(HashProcessor, self).__init__(*args, **kwargs)

    def _compile_results(self):
        pass

    def _handle_new_result(self, result):
        pass


class StatementProcessor(Processor):
    """A Processor to get Statements from the server."""
    def __init__(self, query: Query, limit=None, sort_by='ev_count', ev_limit=10,
                 filter_ev=True, timeout=None, strict_stop=False, persist=True,
                 use_obtained_counts=False):

        self.statements = []
        self.statements_sample = None

        self.__statement_jsons = {}
        self.__evidence_counts = {}
        self.__source_counts = {}

        self.use_obtained_counts = use_obtained_counts
        self._set_special_params(ev_limit=ev_limit, filter_ev=filter_ev)
        super(StatementProcessor, self).\
            __init__(query, limit=limit, sort_by=sort_by, timeout=timeout,
                     strict_stop=strict_stop, persist=persist)

    def _handle_new_result(self, result):
        """Merge these statement jsons with new jsons."""
        # Merge counts.
        self.__evidence_counts.update(standardize_counts(result.evidence_counts))
        self.__source_counts.update(standardize_counts(result.source_counts))

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


class AgentProcessor(Processor):
    """A Processor to get Agent pairs from the server."""
    def __init__(self, query: Query, limit=None, sort_by='ev_count',
                 with_hashes=False, timeout=None, strict_stop=False,
                 persist=True):
        self.complexes_covered = set()
        self._set_special_params(with_hashes=with_hashes,
                                 complexes_covered=self.complexes_covered)
        super(AgentProcessor, self).\
            __init__(query, limit=limit, sort_by=sort_by, timeout=timeout,
                     strict_stop=strict_stop, persist=persist)

    def _compile_results(self):
        pass

    def _handle_new_result(self, result):
        pass
