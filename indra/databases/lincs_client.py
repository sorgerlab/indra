from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_drug_target_data', 'get_small_molecule_data',
           'get_protein_data', 'LincsClient']

import csv
import requests
from os import path

LINCS_URL = 'http://lincs.hms.harvard.edu/db'


class LincsClient(object):
    def __init__(self):
        self._sm_data = get_small_molecule_data()
        self._prot_data = get_protein_data()
        return

    def get_small_molecule_name(self, id_val, id_type='hms-lincs'):
        """Get the name of a small molecule from the LINCS sm metadata.

        Parameters
        ----------
        id_val : str
            The value of the id
        id_type : str
            An indication of the type of id. Handled options are: 'hms-lincs',
            'short-hms-lincs', and 'chembl'.

        Returns
        -------
        str
            The name of the small molecule.
        """
        res = self.__harvest_sm_data(lambda _, info_dict: info_dict['Name'],
                                     id_val, id_type)
        # To end up with a name, here we take the value of the result dict.
        # Note that in every case in this dataset, there is only
        # one unique value in the dict i.e. len(set(res.values())) == 1
        # so it is safe to do this.
        if isinstance(res, dict):
            return list(res.values())[0]
        else:
            return res


    def get_small_molecule_ref(self, id_val, id_type='hms-lincs'):
        """Get the id refs of a small molecule from the LINCS sm metadata.

        Parameters
        ----------
        id_val : str
            The value of the id
        id_type : str
            An indication of the type of id. Handles options are: 'hms-lincs',
            'short-hms-lincs', and 'chembl'.

        Returns
        -------
        The returns vary. If the id is hms-lincs, it will return a single value,
        however if there are multiple entries for the value, it will return a
        dict of values keyed by hms-lincs ids.
        """
        mappings = dict(chembl='ChEMBL ID', chebi='ChEBI ID',
                        pubchem='PubChem CID', inchi='InChi Key',
                        lincs='LINCS ID')

        def ret_func(lincs_id, info_dict):
            return _build_db_refs(lincs_id, info_dict, **mappings)

        return self.__harvest_sm_data(ret_func, id_val, id_type)

    def __harvest_sm_data(self, ret_func, id_val, id_type):
        if id_type not in ['hms-lincs', 'short-hms-lincs', 'chembl']:
            raise ValueError('Unexpected value for input id_type: %s'
                             % id_type)

        if id_type == 'short-hms-lincs':
            ret = {}
            for hms_lincs_id, info_dict in self._sm_data.items():
                if hms_lincs_id.split('-')[0] == id_val:
                    ret[hms_lincs_id] = ret_func(hms_lincs_id, info_dict)
        elif id_type == 'chembl':
            ret = {}
            for hms_lincs_id, info_dict in self._sm_data.items():
                if info_dict['ChEMBL ID'] == id_val:
                    ret[hms_lincs_id] = ret_func(hms_lincs_id, info_dict)
        else:
            ret = ret_func(id_val, self._sm_data[id_val])
        return ret

    def get_protein_ref(self, id_val, id_type='hms-lincs'):
        """Get the refs for a protein from the LINCs protein metadata.

        Parameters
        ----------
        id_val : str
            The value of the id
        id_type : str
            An indication of the type of id. Handles options are: 'hms-lincs'
            and 'entrez'.
        """
        if id_type not in ['hms-lincs', 'entrez']:
            raise ValueError("Unexpected value for input id_type: %s" % id_type)

        # TODO: We could get phosphorylation states from the prtein data.
        if id_type == 'hms-lincs':
            return _build_db_refs(id_val, self._prot_data[id_val],
                                  entrez='Gene ID', uniprot='UniProt ID')
        elif id_type == 'entrez':
            ret = {}
            for hms_id, info_dict in self._prot_data.items():
                if info_dict['Gene ID'] != id_val:
                    continue
                ret[hms_id] = _build_db_refs(hms_id, info_dict,
                                             entrez='Gene ID',
                                             uniprot='UniProt ID')
            return ret


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
