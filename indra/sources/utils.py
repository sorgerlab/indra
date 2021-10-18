# -*- coding: utf-8 -*-

"""Processor for remote INDRA JSON files."""

import pickle
from typing import ClassVar, Iterable, List, TYPE_CHECKING

import requests

from ..statements import Statement, print_stmt_summary, stmts_from_json

if TYPE_CHECKING:
    import click

__all__ = [
    "Processor",
    'RemoteProcessor',
    "SimpleProcessor",
]


class Processor:
    """A base class for processors."""

    name: ClassVar[str]

    def extract_statements(self) -> List[Statement]:
        """Extract statements from the remote JSON file."""
        raise NotImplementedError

    @classmethod
    def get_cli(cls) -> "click.Command":
        """Get the CLI for this processor."""
        import click

        @click.command()
        @click.option('--save', is_flag=True)
        def _main(save: bool):
            click.secho(cls.name, fg='green', bold=True)
            inst = cls()
            stmts = inst.extract_statements()
            if save:
                import pystow
                with pystow.join("indra", cls.name, name="stmts.pkl").open("wb") as file:
                    pickle.dump(stmts, file, protocol=pickle.HIGHEST_PROTOCOL)
            print_stmt_summary(stmts)

        return _main

    @classmethod
    def cli(cls) -> None:
        """Run the CLI for this processor."""
        _main = cls.get_cli()
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


class SimpleProcessor(Processor):
    """A processor that is implemented by creating an iterator over statements."""

    #: A list of INDRA statements extracted by the processor
    statements: List[Statement]

    def __init__(self):
        self.statements = []

    def iter_statements(self) -> Iterable[Statement]:
        """Generate statements with an iterable."""
        raise NotImplementedError

    def extract_statements(self) -> List[Statement]:
        """Generate and store statements if not pre-cached, then return then."""
        if not self.statements:
            self.statements = list(self.iter_statements())
        return self.statements
