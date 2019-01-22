from __future__ import absolute_import, unicode_literals
from builtins import dict, str

from indra.sources.indra_db_rest.util import submit_query_request
from indra.sources.indra_db_rest.exceptions import IndraDBRestResponseError

__all__ = ['IndraDBRestProcessor']

import logging
from threading import Thread
from datetime import datetime
from collections import OrderedDict, defaultdict

from indra.statements import stmts_from_json


logger = logging.getLogger(__name__)


class IndraDBRestProcessor(object):
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
        return self.get_ev_count_by_hash(stmt.get_hash(shallow=True))

    def get_ev_count_by_hash(self, stmt_hash):
        """Get the total evidence count for a statement hash."""
        return self.__evidence_counts.get(str(stmt_hash))

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
        resp = submit_query_request('from_agents', *agent_strs, **params)
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
        """Slightly lower level function gets statements from the REST API."""
        # Handle the content if we were limited.
        args = [agent_strs, stmt_types, params, persist]
        logger.info("The remainder of the query will be performed in a "
                    "thread...")
        self.__th = Thread(target=self._run_queries, args=args)
        self.__th.start()

        if block_secs is None:
            logger.info("Waiting for thread to complete...")
            self.__th.join()
        elif block_secs:  # is not 0
            logger.info("Waiting at most %d seconds for thread to complete..."
                        % block_secs)
            self.__th.join(block_secs)
        return


