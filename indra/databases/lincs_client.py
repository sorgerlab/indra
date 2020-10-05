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
from indra.databases.identifiers import ensure_chembl_prefix


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
        extra_sm_data = load_lincs_extras()
        self._sm_data.update(extra_sm_data)

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
                key = k.upper()
                value = entry[v]
                # Swap in primary PubChem IDs where there is an outdated one
                if key == 'PUBCHEM' and value in pc_to_primary_mappings:
                    value = pc_to_primary_mappings[value]
                # Fix CHEMBL IDs
                if key == 'CHEMBL':
                    value = ensure_chembl_prefix(value)
                refs[key] = value
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


def load_lincs_extras():
    fname = os.path.join(resources, 'hms_lincs_extra.tsv')
    with open(fname, 'r') as fh:
        rows = [line.strip('\n').split('\t') for line in fh.readlines()]
    return {r[0]: {'HMS LINCS ID': r[0],
                   'Name': r[1],
                   'ChEMBL ID': r[2] if r[2] else ''}
            for r in rows[1:]}


# This is a set of mappings specific to HMS-LINCS that map outdated compound
# IDs appearing in HMS-LINCS to preferred compound IDs. This can be obtained
# more generally via indra.databases.pubchem_client, but this is a pre-compiled
# version here for fast lookups in this client.
pc_to_primary_mappings = \
    {'23624255': '135564985',
     '10451420': '135465539',
     '10196499': '135398501',
     '57899889': '135564632',
     '53239990': '135564599',
     '71433937': '136240579',
     '53401173': '135539077',
     '71543332': '135398499',
     '5353940': '5169',
     '49830557': '135398510',
     '11258443': '135451019',
     '68925359': '135440466',
     '16750408': '135565545',
     '57347681': '135565635',
     '5357795': '92577',
     '56965966': '135398516',
     '24906282': '448949',
     '66524294': '135398492',
     '11696609': '135398495',
     '9549301': '135473382',
     '56965894': '135423438',
     }
