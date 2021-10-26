# -*- coding: utf-8 -*-

"""API for CREEDS."""

import json
from pathlib import Path
from typing import Union

from .processor import (
    CREEDSChemicalProcessor,
    CREEDSDiseaseProcessor,
    CREEDSGeneProcessor,
    CREEDSProcessor,
)

processors = {
    "gene": CREEDSGeneProcessor,
    "disease": CREEDSDiseaseProcessor,
    "chemical": CREEDSChemicalProcessor,
}

__all__ = [
    "process_from_file",
    "process_from_web",
]


def process_from_web(entity_type: str) -> CREEDSProcessor:
    """Process statements from CREEDS by automatially downloading them.

    Parameters
    ----------
    entity_type :
        Either 'gene', 'disease', or 'chemical' to specify
        which dataset to get.

    Returns
    -------
    processor : SimpleProcessor
        A processor with pre-extracted statements.
    """
    return _load(entity_type)


def process_from_file(
    path: Union[str, Path],
    entity_type: str,
) -> CREEDSProcessor:
    """Process statements from CREEDS in a file.

    Parameters
    ----------
    path :
        The path to a JSON file containing records for the CREEDS data
    entity_type :
        Either 'gene', 'disease', or 'chemical' to specify
        which dataset to get.

    Returns
    -------
    processor : SimpleProcessor
        A processor with pre-extracted statements.
    """
    with open(path) as file:
        records = json.load(file)
    if len(records) != 1:
        raise ValueError
    return _load(entity_type, records)


def _load(entity_type, records=None):
    processor_cls = processors[entity_type]
    processor = processor_cls(records)
    processor.extract_statements()
    return processor
