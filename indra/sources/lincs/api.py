from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = []

import requests

DATASET_URL = 'http://lincs.hms.harvard.edu/db/datasets/20000/results'


def _get_lincs_drug_target_data():
    resp = requests.get(DATASET_URL, params={'output_type': '.csv'})
    assert resp.status_code == 200, resp.text
    csv_str = resp.content.decode('utf-8')
    return [tuple(line.split(',')) for line in csv_str.splitlines()]
