import json
import os
import re
from functools import lru_cache
from urllib.parse import urlencode
from os.path import abspath, dirname, join, pardir
import requests
from indra.util import read_unicode_csv

MESH_URL = 'https://id.nlm.nih.gov/mesh/'
HERE = dirname(abspath(__file__))
RESOURCES = join(HERE, pardir, 'resources')
MESH_FILE = join(RESOURCES, 'mesh_id_label_mappings.tsv')
MESH_SUPP_FILE = join(RESOURCES, 'mesh_supp_id_label_mappings.tsv')
GO_MAPPINGS = join(RESOURCES, 'mesh_go_mappings.tsv')


mesh_id_to_name = {}
mesh_name_to_id = {}
mesh_name_to_id_name = {}


def _load_mesh_file(path):
    it = read_unicode_csv(path, delimiter='\t')
    for mesh_id, mesh_label, mesh_terms_str in it:
        mesh_id_to_name[mesh_id] = mesh_label
        mesh_name_to_id[mesh_label] = mesh_id
        mesh_terms = mesh_terms_str.split('|')
        for term in mesh_terms:
            mesh_name_to_id_name[term] = [mesh_id, mesh_label]


_load_mesh_file(MESH_FILE)
if os.path.exists(MESH_SUPP_FILE):
    _load_mesh_file(MESH_SUPP_FILE)


def _load_go_mappings(path):
    mesh_to_go = {}
    go_to_mesh = {}
    for mesh_id, go_id in read_unicode_csv(path, delimiter='\t'):
        mesh_to_go[mesh_id] = go_id
        go_to_mesh[go_id] = mesh_id
    return mesh_to_go, go_to_mesh


mesh_to_go, go_to_mesh = _load_go_mappings(GO_MAPPINGS)


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
    url = MESH_URL + mesh_id + '.json'
    resp = requests.get(url)
    if resp.status_code != 200:
        return None
    mesh_json = resp.json()
    try:
        label = mesh_json['label']['@value']
    except (KeyError, IndexError, TypeError) as e:
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
        Whether to allow queries to the NLM REST API if the given MESH ID is
        not contained in INDRA's internal MESH mappings file. Default is False
        (allows REST API queries).

    Returns
    -------
    str
        Label for the MESH ID, or None if the query failed or no label was
        found.
    """
    indra_mesh_mapping = mesh_id_to_name.get(mesh_id)
    if offline or indra_mesh_mapping is not None:
        return indra_mesh_mapping
    # Look up the MESH mapping from NLM if we don't have it locally
    return get_mesh_name_from_web(mesh_id)


def get_mesh_id_name(mesh_term, offline=False):
    """Get the MESH ID and name for the given MESH term.

    Uses the mappings table in `indra/resources`; if the MESH term is not
    listed there, falls back on the NLM REST API.

    Parameters
    ----------
    mesh_term : str
        MESH Descriptor or Concept name, e.g. 'Breast Cancer'.
    offline : bool
        Whether to allow queries to the NLM REST API if the given MESH term is
        not contained in INDRA's internal MESH mappings file. Default is False
        (allows REST API queries).

    Returns
    -------
    tuple of strs
        Returns a 2-tuple of the form `(id, name)` with the ID of the
        descriptor corresponding to the MESH label, and the descriptor name
        (which may not exactly match the name provided as an argument if it is
        a Concept name). If the query failed, or no descriptor corresponding to
        the name was found, returns a tuple of (None, None).
    """
    indra_mesh_id = mesh_name_to_id.get(mesh_term)
    if indra_mesh_id is not None:
        return indra_mesh_id, mesh_term

    indra_mesh_id, new_term = \
        mesh_name_to_id_name.get(mesh_term, (None, None))
    if indra_mesh_id is not None:
        return indra_mesh_id, new_term

    if offline:
        return None, None

    # Look up the MESH mapping from NLM if we don't have it locally
    return get_mesh_id_name_from_web(mesh_term)


@lru_cache(maxsize=1000)
def submt_sparql_query(query_body):
    url = MESH_URL + 'sparql'
    query = '%s\n%s' % (mesh_rdf_prefixes, query_body)
    args = {'query': query, 'format': 'JSON', 'inference': 'true'}
    resp = requests.get(url, params=args)
    # Check status
    if resp.status_code != 200:
        return None
    try:
        # Try to parse the json response (this can raise exceptions if we
        # got no response).
        return resp.json()
    except Exception:
        return None


def get_mesh_id_name_from_web(mesh_term):
    """Get the MESH ID and name for the given MESH term using the NLM REST API.

    Parameters
    ----------
    mesh_term : str
        MESH Descriptor or Concept name, e.g. 'Breast Cancer'.

    Returns
    -------
    tuple of strs
        Returns a 2-tuple of the form `(id, name)` with the ID of the
        descriptor corresponding to the MESH label, and the descriptor name
        (which may not exactly match the name provided as an argument if it is
        a Concept name). If the query failed, or no descriptor corresponding to
        the name was found, returns a tuple of (None, None).
    """
    query_body = """
        SELECT ?d ?dName ?c ?cName
        FROM <http://id.nlm.nih.gov/mesh>
        WHERE {
          ?d a meshv:Descriptor .
          ?d meshv:concept ?c .
          ?d rdfs:label ?dName .
          ?c rdfs:label ?cName
          FILTER (REGEX(?dName,'^%s$','i') || REGEX(?cName,'^%s$','i'))
        }
        ORDER BY ?d
    """ % (mesh_term, mesh_term)
    mesh_json = submt_sparql_query(query_body)
    if mesh_json is None:
        return None, None
    try:
        # Choose the first entry (should usually be only one)
        id_uri = mesh_json['results']['bindings'][0]['d']['value']
        name = mesh_json['results']['bindings'][0]['dName']['value']
    except (KeyError, IndexError, json.decoder.JSONDecodeError) as e:
        return None, None

    # Strip the MESH prefix off the ID URI
    m = re.match('http://id.nlm.nih.gov/mesh/([A-Za-z0-9]*)', id_uri)
    assert m is not None
    id = m.groups()[0]
    return id, name


def mesh_isa(mesh_id1, mesh_id2):
    query_body = """
        SELECT DISTINCT ?o
        FROM <http://id.nlm.nih.gov/mesh>
        WHERE {
          mesh:%s meshv:broaderDescriptor+ ?o .
        }
        """ % mesh_id1
    mesh_json = submt_sparql_query(query_body)
    if mesh_json is None:
        return False
    try:
        results = mesh_json['results']['bindings']
        for result in results:
            id_uri = result['o']['value']
            # Strip the MESH prefix off the ID URI
            m = re.match('http://id.nlm.nih.gov/mesh/([A-Za-z0-9]*)', id_uri)
            id = m.groups()[0]
            if mesh_id2 == id:
                return True
        return False
    except Exception:
        return False


def get_mesh_tree_number(mesh_id):
    query_body = """
        SELECT DISTINCT ?tn
        FROM <http://id.nlm.nih.gov/mesh>
        WHERE {
          mesh:%s meshv:treeNumber ?tn
        }
        """ % mesh_id
    mesh_json = submt_sparql_query(query_body)
    if mesh_json is None:
        return None
    try:
        results = mesh_json['results']['bindings']
        tree_uri = results[0]['tn']['value']
        m = re.match('http://id.nlm.nih.gov/mesh/([A-Z0-9.]*)', tree_uri)
        tree = m.groups()[0]
        return tree
    except Exception:
        return None


def get_go_id(mesh_id):
    """Return a GO ID corresponding to the given MeSH ID.

    Parameters
    ----------
    mesh_id : str
        MeSH ID to map to GO

    Returns
    -------
    str
        The GO ID corresponding to the given MeSH ID, or None if not available.
    """
    return mesh_to_go.get(mesh_id)


def get_mesh_id_from_go_id(go_id):
    """Return a MeSH ID corresponding to the given GO ID.

    Parameters
    ----------
    go_id : str
        GO ID to map to MeSH

    Returns
    -------
    str
        The MeSH ID corresponding to the given GO ID, or None if not
        available.
    """
    return go_to_mesh.get(go_id)


mesh_rdf_prefixes = """
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
        PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
        PREFIX mesh2019: <http://id.nlm.nih.gov/mesh/2019/>
        PREFIX mesh2018: <http://id.nlm.nih.gov/mesh/2018/>
        PREFIX mesh2017: <http://id.nlm.nih.gov/mesh/2017/>
    """
