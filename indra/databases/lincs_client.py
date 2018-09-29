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

    def get_small_molecule_name(self, hms_lincs_id):
        """Get the name of a small molecule from the LINCS sm metadata.

        Parameters
        ----------
        hms_lincs_id : str
            The HMS LINCS ID of the small molecule.

        Returns
        -------
        str
            The name of the small molecule.
        """
        entry = self._get_entry_by_id(hms_lincs_id)
        if not entry:
            return None
        name = entry['Name']
        return name

    def get_small_molecule_refs(self, hms_lincs_id):
        """Get the id refs of a small molecule from the LINCS sm metadata.

        Parameters
        ----------
        hms_lincs_id : str
            The HMS LINCS ID of the small molecule.

        Returns
        -------
        dict
            A dictionary of references.
        """
        refs = {'HMS-LINCS': hms_lincs_id}

        entry = self._get_entry_by_id(hms_lincs_id)
        # If there is no entry for this ID
        if not entry:
            return refs

        # If there is an entry then fill up the refs with existing values
        mappings = dict(chembl='ChEMBL ID', chebi='ChEBI ID',
                        pubchem='PubChem CID', lincs='LINCS ID')
        for k, v in mappings.items():
            if entry.get(v):
                refs[k.upper()] = entry.get(v)
        return refs

    def _get_entry_by_id(self, hms_lincs_id):
        # This means it's a short ID
        if '-' not in hms_lincs_id:
            keys = [k for k in self._sm_data.keys() if
                    k.startswith(hms_lincs_id)]
            if not keys:
                logger.error('Couldn\'t find entry for %s' % hms_lincs_id)
                return None
            entry = self._sm_data[keys[0]]
        # This means it's a full ID
        else:
            entry = self._sm_data.get(hms_lincs_id)
            if not entry:
                logger.error('Couldn\'t find entry for %s' % hms_lincs_id)
                return None
        return entry

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
                                  egid='Gene ID', up='UniProt ID')
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
    print(db_refs)
    return db_refs
