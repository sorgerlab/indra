# -*- coding: utf-8 -*-

"""Processor for remote INDRA JSON files."""

import requests
from typing import List

from ..statements import Statement, stmts_from_json

__all__ = [
    'RemoteProcessor',
]


class RemoteProcessor:
    """A processor for INDRA JSON file to be retrieved by URL."""

    #: The URL of the data
    url: str

    #: A list of statements
    statements: List[Statement]

    def __init__(self, url: str):
        self.url = url
        self.statements = []

    def extract_statements(self) -> List[Statement]:
        """Extract statements from the remote JSON file."""
        res = requests.get(self.url)
        res.raise_for_status()
        self.statements = stmts_from_json(res.json())
        return self.statements
