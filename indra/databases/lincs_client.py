from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_drug_target_data', 'get_small_molecule_data',
           'get_protein_data']

import csv
import requests
from os import path

LINCS_URL = 'http://lincs.hms.harvard.edu/db'


class LincsClient(object):
    def __init__(self):
        self._sm_data = get_small_molecule_data()
        self._prot_data = get_protein_data()
        return

    def get_small_molecule_ref(self, lincs_id):
        return _build_db_refs(lincs_id, self._sm_data[lincs_id],
                              chembl='ChEMBL ID', chebi='ChEBI ID',
                              pubchem='PubChem CID', inchi='InChi Key',
                              lincs='LINCS ID')

    def get_protein_ref(self, lincs_id):
        # TODO: We could get phosphorylation states from the prtein data.
        return _build_db_refs(lincs_id, self._prot_data[lincs_id],
                              hgnc='Gene ID', uniprot='UniProt ID')


def get_drug_target_data():
    """Load the csv into a list of dicts containing the LINCS drug target data.

    Returns
    -------
    data : list[dict]
        A list of dicts, each keyed based on the header of the csv, with values
        as the corresponding column values.
    """
    url = path.join(LINCS_URL, 'datasets/20000/results')
    return _load_lincs_csv(url)


def get_small_molecule_data():
    """Load the csv of LINCS small molecule metadata into a dict.

    Returns
    -------
    sm_dict : dict[dict]
        A dict keyed by HMS LINCS small molecule ids, with the metadata
        contained in a dict of row values keyed by the column headers extracted
        from the csv.
    """
    url = path.join(LINCS_URL, 'sm/')  # The trailing / is deliberate
    sm_data = _load_lincs_csv(url)
    sm_dict = {d['HMS LINCS ID']: d.copy() for d in sm_data}
    assert len(sm_dict) == len(sm_data), "We lost data."
    return sm_dict


def get_protein_data():
    """Load the csv of LINCS protein metadata into a dict.

    Returns
    -------
    prot_dict : dict[dict]
        A dict keyed by HMS LINCS protein ids, with the metadata contained in a
        dict of row values keyed by the column headers extracted from the csv.
    """
    url = path.join(LINCS_URL, 'proteins/')
    prot_data = _load_lincs_csv(url)
    prot_dict = {d['HMS LINCS ID']: d.copy() for d in prot_data}
    assert len(prot_dict) == len(prot_data), "We lost data."
    return prot_dict


def _load_lincs_csv(url):
    """Helper function to turn csv rows into dicts."""
    resp = requests.get(url, params={'output_type': '.csv'})
    assert resp.status_code == 200, resp.text
    csv_str = resp.content.decode('utf-8')
    csv_lines = csv_str.splitlines()
    headers = csv_lines[0].split(',')
    return [{headers[i]: val for i, val in enumerate(line_elements)}
            for line_elements in csv.reader(csv_lines[1:])]


def _build_db_refs(lincs_id, data, **mappings):
    db_refs = {'HMS-LINCS': lincs_id}
    for db_ref, key in mappings.items():
        if data[key]:
            db_refs[db_ref.upper()] = data[key]
    return db_refs
