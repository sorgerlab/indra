from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from collections import namedtuple
from indra.util import read_unicode_csv

logger = logging.getLogger('signor')


_signor_fields = [
    'ENTITYA',
    'TYPEA',
    'IDA',
    'DATABASEA',
    'ENTITYB',
    'TYPEB',
    'IDB',
    'DATABASEB',
    'EFFECT',
    'MECHANISM',
    'RESIDUE',
    'SEQUENCE',
    'TAX_ID',
    'CELL_DATA',
    'TISSUE_DATA',
    'MODULATOR_COMPLEX',
    'TARGET_COMPLEX',
    'MODIFICATIONA',
    'MODASEQ',
    'MODIFICATIONB',
    'MODBSEQ',
    'PMID',
    'DIRECT',
    'NOTES',
    'ANNOTATOR',
    'SENTENCE',
    'SIGNOR_ID',
]


SignorRow = namedtuple('SignorRow', _signor_fields)


class SignorProcessor(object):
    """Processor for Signor dataset available at http://signor.uniroma2.it.

    See publication:

    Perfetto et al., "SIGNOR: a database of causal relationships between
    biological entities," Nucleic Acids Research, Volume 44, Issue D1, 4
    January 2016, Pages D548–D554. https://doi.org/10.1093/nar/gkv1048

    Parameters
    ----------
    signor_csv : str
        Path to SIGNOR CSV file.
    delimiter : str
        Field delimiter for CSV file. Defaults to semicolon ';'.
    """
    def __init__(self, signor_csv, delimiter=';'):
        # Get generator over the CSV file
        data_iter = read_unicode_csv(signor_csv, delimiter=';')
        # Skip the header row
        next(data_iter)
        # Process into a list of SignorRow namedtuples
        self._data = [SignorRow(*r) for r in data_iter]

    def _process_row(self):
        pass


