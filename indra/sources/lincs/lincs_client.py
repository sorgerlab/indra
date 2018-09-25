from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_lincs_drug_target_data']

import requests
from os import path

LINCS_URL = 'http://lincs.hms.harvard.edu/db'


def get_lincs_drug_target_data():
    url = path.join(LINCS_URL, '/datasets/20000/results')
    return _load_lincs_csv(url)


def get_lincs_small_molecule_data():
    url = path.join(LINCS_URL, 'sm/')  # The trailing / is deliberate
    return _load_lincs_csv(url)


def _load_lincs_csv(url):
    resp = requests.get(url, params={'output_type', '.csv'})
    assert resp.status_code == 200, resp.text
    csv_str = resp.content.decode('utf-8')
    csv_lines = csv_str.splitlines()
    headers = csv_lines[0].split(',')
    return [{headers[i]: val for i, val in enumerate(line_elements)}
            for line_elements in csv.reader(csv_lines[1:])]
