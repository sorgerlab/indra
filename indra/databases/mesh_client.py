import requests
import functools

mesh_url = 'http://id.nlm.nih.gov/mesh/'


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

