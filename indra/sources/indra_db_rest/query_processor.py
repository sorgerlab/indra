import json
import requests

from indra.config import get_config
from .query import Query


class Processor:
    def __init__(self, query: Query, limit=None, sort_by='default'):
        self.query = query
        self.limit = limit
        self.sort_by = sort_by
        self._run()

    def _run(self):
        raise NotImplementedError()

    def _make_query(self):
        """Make a single query to the server."""
        query_json = self.query.to_json()
        url = get_config('INDRA_DB_REST_URL', failure_ok=False)
        resp = requests.post(url, data=json.dumps(query_json))


class HashProcessor(Processor):
    """A processor to get hashes from the server."""

    def __init__(self, *args, **kwargs):
        self.hashes = []
        super(HashProcessor, self).__init__(*args, **kwargs)

    def _run(self):
        pass


class StatementProcessor(Processor):
    """A Processor to get Statements from the server."""
    def _run(self):
        pass


class AgentProcessor(Processor):
    """A Processor to get Agent pairs from the server."""
    def _run(self):
        pass
