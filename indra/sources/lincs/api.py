from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['process_from_web']

from indra.sources.lincs.processor import LincsProcessor

from .lincs_client import get_lincs_drug_target_data


def process_from_web():
    lincs_data = get_lincs_drug_target_data()
    return LincsProcessor(lincs_data)
