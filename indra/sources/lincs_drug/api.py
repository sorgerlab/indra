from __future__ import absolute_import, print_function, unicode_literals

__all__ = ['process_from_web']

from .processor import LincsProcessor
from indra.databases.lincs_client import get_drug_target_data


def process_from_web():
    """Return a processor for the LINCS drug target data.

    Returns
    -------
    LincsProcessor
        A LincsProcessor object which contains extracted INDRA Statements
        in its statements attribute.
    """
    lincs_data = get_drug_target_data()
    return LincsProcessor(lincs_data)
