# -*- coding: utf-8 -*-

"""API for CREEDS."""

from .processor import CREEDSChemicalProcessor, CREEDSDiseaseProcessor, CREEDSGeneProcessor
from ..utils import SimpleProcessor

processors = {
    'gene': CREEDSGeneProcessor,
    'disease': CREEDSDiseaseProcessor,
    'chemical': CREEDSChemicalProcessor,
}

__all__ = [
    'process_from_web',
]


def process_from_web(entity_type: str) -> SimpleProcessor:
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
    processor_cls = processors[entity_type]
    processor = processor_cls()
    processor.extract_statements()
    return processor
