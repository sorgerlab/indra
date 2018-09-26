from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['process_from_web']

from .processor import LincsProcessor
from .lincs_client import get_drug_target_data


def process_from_web():
    """Get a processor for the LINCS drug target data."""
    lincs_data = get_drug_target_data()
    return LincsProcessor(lincs_data)
