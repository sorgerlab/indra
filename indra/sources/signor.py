from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.util import read_unicode_csv


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
        self._data = list(read_unicode_csv(signor_csv, delimiter=';'))
        

