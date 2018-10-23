import csv
import requests

mesh_url = 'http://id.nlm.nih.gov/mesh/'

# Python3
try:
    from functools import lru_cache
# Python2
except ImportError:
    from functools32 import lru_cache


def 

def get_mesh_name(mesh_id):
    """Get the MESH label for the given MESH ID.

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


if __name__ == '__main__':
    with open('mesh_ids.txt', 'rt') as f:
        mesh_ids = [line.strip() for line in f.readlines()]
    mesh_mappings = []
    for mid in mesh_ids:
        if mid.startswith('MESH:'):
            mid = mid[5:]
        mesh_label = get_mesh_name(mid)
        mesh_mappings.append((mid, mesh_label))
        print(f"{mid} -> {mesh_label}")
    with open('mesh_id_label_mappings.tsv', 'wt') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        csvwriter.writerows(mesh_mappings)

