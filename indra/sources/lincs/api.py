from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = []

import csv
import requests

from indra.sources.lincs.processor import LincsProcessor

DATASET_URL = 'http://lincs.hms.harvard.edu/db/datasets/20000/results'


def process_from_web():
    lincs_data = _get_lincs_drug_target_data()
    return LincsProcessor(lincs_data)


def _get_lincs_drug_target_data():
    resp = requests.get(DATASET_URL, params={'output_type': '.csv'})
    assert resp.status_code == 200, resp.text
    csv_str = resp.content.decode('utf-8')
    csv_lines = csv_str.splitlines()
    headers = csv_lines[0].split(',')
    return [{headers[i]: val for i, val in enumerate(line_elements)}
            for line_elements in csv.reader(csv_lines[1:])]
