from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_drug_target_data', 'LincsClient', 'load_lincs_csv']

import os
import sys
import json
import logging
import requests
from io import StringIO, BytesIO
from indra.util import read_unicode_csv_fileobj


logger = logging.getLogger(__name__)


LINCS_URL = 'http://lincs.hms.harvard.edu/db'


resources = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                         os.path.pardir, 'resources')
lincs_sm = os.path.join(resources, 'lincs_small_molecules.json')
lincs_prot = os.path.join(resources, 'lincs_proteins.json')


class LincsClient(object):
    """Client for querying LINCS small molecules and proteins."""
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
        entry = self._get_entry_by_id(self._sm_data, hms_lincs_id)
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

        entry = self._get_entry_by_id(self._sm_data, hms_lincs_id)
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

    def get_protein_refs(self, hms_lincs_id):
        """Get the refs for a protein from the LINCs protein metadata.

        Parameters
        ----------
        hms_lincs_id : str
            The HMS LINCS ID for the protein

        Returns
        -------
        dict
            A dictionary of protein references.
        """
        # TODO: We could get phosphorylation states from the protein data.
        refs = {'HMS-LINCS': hms_lincs_id}

        entry = self._get_entry_by_id(self._prot_data, hms_lincs_id)
        # If there is no entry for this ID
        if not entry:
            return refs
        mappings = dict(egid='Gene ID', up='UniProt ID')
        for k, v in mappings.items():
            if entry.get(v):
                refs[k.upper()] = entry.get(v)
        return refs

    def _get_entry_by_id(self, resource, hms_lincs_id):
        # This means it's a short ID
        if '-' not in hms_lincs_id:
            keys = [k for k in resource.keys() if k.startswith(hms_lincs_id)]
            if not keys:
                logger.debug('Couldn\'t find entry for %s' % hms_lincs_id)
                return None
            entry = resource[keys[0]]
        # This means it's a full ID
        else:
            entry = resource.get(hms_lincs_id)
            if not entry:
                logger.debug('Couldn\'t find entry for %s' % hms_lincs_id)
                return None
        return entry


def get_drug_target_data():
    """Load the csv into a list of dicts containing the LINCS drug target data.

    Returns
    -------
    data : list[dict]
        A list of dicts, each keyed based on the header of the csv, with values
        as the corresponding column values.
    """
    url = LINCS_URL + '/datasets/20000/results'
    return load_lincs_csv(url)


def _build_db_refs(lincs_id, data, **mappings):
    db_refs = {'HMS-LINCS': lincs_id}
    for db_ref, key in mappings.items():
        if data[key]:
            db_refs[db_ref.upper()] = data[key]
    return db_refs


def load_lincs_csv(url):
    """Helper function to turn csv rows into dicts."""
    resp = requests.get(url, params={'output_type': '.csv'}, timeout=120)
    resp.raise_for_status()
    if sys.version_info[0] < 3:
        csv_io = BytesIO(resp.content)
    else:
        csv_io = StringIO(resp.text)
    data_rows = list(read_unicode_csv_fileobj(csv_io, delimiter=','))
    headers = data_rows[0]
    return [{header: val for header, val in zip(headers, line_elements)}
            for line_elements in data_rows[1:]]

