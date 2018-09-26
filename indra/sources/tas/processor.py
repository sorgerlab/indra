from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = []

import csv
from os import path

HERE = path.dirname(path.abspath(__file__))
DATAFILE_NAME = 'classification_hms_cmpds_symbol.csv'


def _load_data():
    csv_path = path.join(HERE, path.pardir, path.pardir, path.pardir, 'data',
                         DATAFILE_NAME)
    with open(csv_path, 'r') as f:
        csv_lines = f.readlines()
    reader = csv.reader(csv_lines)
    headers = reader.__next__()
    return [{headers[i]: val for i, val in enumerate(line)} for line in reader]
