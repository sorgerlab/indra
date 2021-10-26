# -*- coding: utf-8 -*-

"""API for CREEDS."""

import json
import requests
from pathlib import Path
from typing import Union

from .processor import (
    CREEDSChemicalProcessor,
    CREEDSDiseaseProcessor,
    CREEDSGeneProcessor,
    CREEDSProcessor,
)

BASE_URL = "http://amp.pharm.mssm.edu/CREEDS/download"
urls = {
    "gene": f"{BASE_URL}/single_gene_perturbations-v1.0.json",
    "disease": f"{BASE_URL}/disease_signatures-v1.0.json",
    "chemical": f"{BASE_URL}/single_drug_perturbations-v1.0.json",
}

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
    :
        A processor with pre-extracted statements.
    """
    url = urls[entity_type]
    res = requests.get(url)
    res.raise_for_status()
    records = res.json()
    return process_records(records, entity_type)


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
    :
        A processor with pre-extracted statements.
    """
    with open(path) as file:
        records = json.load(file)
    if len(records) != 1:
        raise ValueError
    return process_records(records, entity_type)


def process_records(records, entity_type):
    """Process statements from CREEDS records.

    Parameters
    ----------
    records :
        A list of records from the CREEDS data
    entity_type :
        Either 'gene', 'disease', or 'chemical' to specify
        which dataset the records represent.

    Returns
    -------
    :
        A processor with pre-extracted statements.
    """
    processor_cls = processors[entity_type]
    processor = processor_cls(records)
    processor.extract_statements()
    return processor
