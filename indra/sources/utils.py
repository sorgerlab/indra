# -*- coding: utf-8 -*-

"""Processor for remote INDRA JSON files."""

from collections import Counter

import requests
from typing import Collection, List

from ..statements import Statement, stmts_from_json

__all__ = [
    'print_statement_summary',
    'RemoteProcessor',
]


def print_statement_summary(statements: Collection[Statement]):
    """Print a summary of the statements."""
    from tabulate import tabulate
    print(tabulate(
        Counter(
            statement.__class__.__name__
            for statement in statements
        ).most_common(),
        headers=["Statement Type", "Count"],
        tablefmt='github',
    ))


class RemoteProcessor:
    """A processor for INDRA JSON file to be retrieved by URL.

    Parameters
    ----------
    url :
        The URL of the INDRA JSON file to load
    """

    #: The URL of the data
    url: str

    def __init__(self, url: str):
        self.url = url
        self._statements = None

    @property
    def statements(self) -> List[Statement]:
        """The extracted statements."""
        if self._statements is None:
            self.extract_statements()
        return self._statements

    def extract_statements(self) -> List[Statement]:
        """Extract statements from the remote JSON file."""
        res = requests.get(self.url)
        res.raise_for_status()
        self._statements = stmts_from_json(res.json())
        return self._statements

    def print_summary(self) -> None:
        """Print a summary of the statements."""
        print_statement_summary(self.statements)
