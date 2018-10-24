from os.path import abspath, dirname, join
import csv
import requests
from indra.util import read_unicode_csv

mesh_url = 'http://id.nlm.nih.gov/mesh/'
mesh_file = join(dirname(abspath(__file__)), '..', 'resources',
                 'mesh_id_label_mappings.tsv')

# Python3
try:
    from functools import lru_cache
# Python2
except ImportError:
    from functools32 import lru_cache

mesh_mappings = {}
for mesh_id, mesh_label in read_unicode_csv(mesh_file, delimiter='\t'):
    mesh_mappings[mesh_id] = mesh_label


@lru_cache(maxsize=1000)
def get_mesh_name_from_web(mesh_id):
    """Get the MESH label for the given MESH ID using the NLM REST API.

    Parameters
    ----------
    mesh_id : str
        MESH Identifier, e.g. 'D003094'.

    Returns
    -------
    str
        Label for the MESH ID, or None if the query failed or no label was
        found.
    """
    url = mesh_url + mesh_id + '.json'
    resp = requests.get(url)
    if resp.status_code != 200:
        return None
    mesh_json = resp.json()
    try:
        label = mesh_json['@graph'][0]['label']['@value']
    except (KeyError, IndexError) as e:
        return None
    return label


def get_mesh_name(mesh_id, offline=False):
    """Get the MESH label for the given MESH ID.

    Uses the mappings table in `indra/resources`; if the MESH ID is not listed
    there, falls back on the NLM REST API.

    Parameters
    ----------
    mesh_id : str
        MESH Identifier, e.g. 'D003094'.
    offline : bool
        Whether to allow queries to the NLM REST API if the given MESH ID is not
        contained in INDRA's internal MESH mappings file. Default is False
        (allows REST API queries).

    Returns
    -------
    str
        Label for the MESH ID, or None if the query failed or no label was
        found.
    """
    indra_mesh_mapping = mesh_mappings.get(mesh_id)
    if offline or indra_mesh_mapping is not None:
        return indra_mesh_mapping
    # Look up the MESH mapping from NLM if we don't have it locally
    return get_mesh_name_from_web(mesh_id)
