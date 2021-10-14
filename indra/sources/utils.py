# -*- coding: utf-8 -*-

"""Processor for remote INDRA JSON files."""

from typing import List

import requests

from ..statements import Statement, print_stmt_summary, stmts_from_json

__all__ = [
    "Processor",
    'RemoteProcessor',
]


class Processor:
    """A base class for processors."""

    def extract_statements(self) -> List[Statement]:
        """Extract statements from the remote JSON file."""
        raise NotImplementedError

    @classmethod
    def cli(cls):
        import click
        @click.command()
        def _main():
            inst = cls()
            stmts = inst.extract_statements()
            print_stmt_summary(stmts)

        _main()


class RemoteProcessor(Processor):
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
        print_stmt_summary(self.statements)
