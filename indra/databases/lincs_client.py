from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_drug_target_data', 'LincsClient']

import os
import csv
import json
import requests


LINCS_URL = 'http://lincs.hms.harvard.edu/db'


resources = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                         os.path.pardir, 'resources')
lincs_sm = os.path.join(resources, 'lincs_small_molecules.json')
lincs_prot = os.path.join(resources, 'lincs_proteins.json')


class LincsClient(object):
    def __init__(self):
        with open(lincs_sm, 'r') as fh:
            self._sm_data = json.load(fh)
        with open(lincs_prot, 'r') as fh:
            self._prot_data = json.load(fh)

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

        # TODO: We could get phosphorylation states from the protein data.
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
    url = LINCS_URL + '/datasets/20000/results'
    return _load_lincs_csv(url)


def _build_db_refs(lincs_id, data, **mappings):
    db_refs = {'HMS-LINCS': lincs_id}
    for db_ref, key in mappings.items():
        if data[key]:
            db_refs[db_ref.upper()] = data[key]
    return db_refs
