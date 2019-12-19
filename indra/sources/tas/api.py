from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['process_csv']

from os import path
from .processor import TasProcessor
from indra.util import read_unicode_csv


HERE = path.dirname(path.abspath(__file__))
DATAFILE_NAME = 'indra_tas.csv'


def _load_data():
    """Load the data from the csv in data.

    The "gene_id" is the Entrez gene id, and the "approved_symbol" is the
    standard gene symbol.

    Returns
    -------
    data : list[dict]
        A list of dicts of row values keyed by the column headers extracted from
        the csv file, described above.
    """
    # Get the cwv reader object.
    csv_path = path.join(HERE, path.pardir, path.pardir, 'resources',
                         DATAFILE_NAME)
    data_iter = list(read_unicode_csv(csv_path))

    # Get the headers.
    headers = data_iter[0]

    # For some reason this heading is oddly formatted and inconsistent with the
    # rest, or with the usual key-style for dicts.
    return [{header: val for header, val in zip(headers, line)}
            for line in data_iter[1:]]


def process_csv(affinity_class_limit=2):
    """Return a TasProcessor for the contents of the csv contained in data.

    Interactions are classified into the following classes based on affinity:
      | 1  -- Kd < 100nM
      | 2  -- 100nM < Kd < 1uM
      | 3  -- 1uM < Kd < 10uM
      | 10 -- Kd > 10uM
    By default, only classes 1 and 2 are extracted but the affinity_class_limit
    parameter can be used to change the upper limit of extracted classes.

    Parameters
    ----------
    affinity_class_limit : Optional[int]
        Defines the highest class of binding affinity that is included in the
        extractions. Default: 2

    Returns
    -------
    TasProcessor
        A TasProcessor object which has a list of INDRA Statements extracted
        from the CSV file representing drug-target inhibitions in its
        statements attribute.
    """
    return TasProcessor(_load_data(), affinity_class_limit)
