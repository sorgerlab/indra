from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['process_csv']

import csv
from os import path
from .processor import TasProcessor

HERE = path.dirname(path.abspath(__file__))
DATAFILE_NAME = 'classification_hms_cmpds_symbol.csv'


def _load_data():
    """Load the data from the csv in data.

    The "gene_id" is the Entrez gene id, and the "approved_symbol" is the
    standard gene symbol. The "hms_id" is the LINCS ID for the drug.

    Returns
    -------
    data : list[dict]
        A list of dicts of row values keyed by the column headers extracted from
        the csv file, described above.
    """
    # Get the cwv reader object.
    csv_path = path.join(HERE, path.pardir, path.pardir, path.pardir, 'data',
                         DATAFILE_NAME)
    with open(csv_path, 'r') as f:
        csv_lines = f.readlines()
    reader = csv.reader(csv_lines)

    # Get the headers.
    headers = reader.__next__()

    # For some reason this heading is oddly formatted and inconsistent with the
    # rest, or with the usual key-style for dicts.
    headers[headers.index('Approved.Symbol')] = 'approved_symbol'
    return [{headers[i]: val for i, val in enumerate(line)} for line in reader]


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
