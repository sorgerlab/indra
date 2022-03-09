"""Utilities for INDRA assemblers."""

from typing import Generic, TypeVar
from abc import ABC, abstractmethod

__all__ = [
  "Assembler",
]

X = TypeVar("X")

class Assembler(ABC, Generic[X]):
    """A base class for INDRA assemblers."""

    @classmethod
    def model_from_stmts(cls, statements: Optional[Statement] = None, **kwargs) -> X:
        assembler = cls(statements, **kwargs)
        return assembler.make_model()

    @abstractmethod
    def make_model(*args, **kwargs) -> X:
        raise NotImplementedError
