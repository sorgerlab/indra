import os
import re
import json
import requests
import itertools
from functools import lru_cache
from os.path import abspath, dirname, join, pardir
from indra.util import read_unicode_csv

MESH_URL = 'https://id.nlm.nih.gov/mesh/'
HERE = dirname(abspath(__file__))
RESOURCES = join(HERE, pardir, 'resources')
MESH_FILE = join(RESOURCES, 'mesh_id_label_mappings.tsv')
MESH_SUPP_FILE = join(RESOURCES, 'mesh_supp_id_label_mappings.tsv')
DB_MAPPINGS = join(RESOURCES, 'mesh_mappings.tsv')


mesh_id_to_name = {}
mesh_name_to_id = {}
mesh_name_to_id_name = {}
mesh_id_to_tree_numbers = {}


def _load_mesh_file(path):
    it = read_unicode_csv(path, delimiter='\t')
    for terms in it:
        if len(terms) == 3:
            mesh_id, mesh_label, mesh_terms_str = terms
        else:
            mesh_id, mesh_label, mesh_terms_str, tree_number_str = terms
            # This is a rare corner case where an entry is outside the
            # tree structure, e.g., D005260, D008297
            if not tree_number_str:
                continue
            mesh_id_to_tree_numbers[mesh_id] = tree_number_str.split('|')
        mesh_terms = mesh_terms_str.split('|') if mesh_terms_str else []
        mesh_id_to_name[mesh_id] = mesh_label
        mesh_name_to_id[mesh_label] = mesh_id
        for term in mesh_terms:
            mesh_name_to_id_name[term] = [mesh_id, mesh_label]


_load_mesh_file(MESH_FILE)
if os.path.exists(MESH_SUPP_FILE):
    _load_mesh_file(MESH_SUPP_FILE)


def _load_db_mappings(path):
    mesh_to_db = {}
    db_to_mesh = {}
    to_db_ambigs = set()
    db_to_ambigs = set()
    for _, mesh_id, _, db_ns, db_id, _ in \
            read_unicode_csv(path, delimiter='\t'):
        # Make sure we don't add any one-to-many mappings
        if mesh_id in mesh_to_db:
            to_db_ambigs.add(mesh_id)
            mesh_to_db.pop(mesh_id, None)
        elif mesh_id not in to_db_ambigs:
            mesh_to_db[mesh_id] = (db_ns, db_id)
        # Make sure we don't add any one-to-many reverse mappings
        if (db_ns, db_id) in db_to_mesh:
            db_to_ambigs.add((db_ns, db_id))
            db_to_mesh.pop((db_ns, db_id), None)
        elif (db_ns, db_id) not in db_to_ambigs:
            db_to_mesh[(db_ns, db_id)] = mesh_id
    return mesh_to_db, db_to_mesh


mesh_to_db, db_to_mesh = _load_db_mappings(DB_MAPPINGS)


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
    if not mesh_term:
        return None, None

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
def submit_sparql_query(query_body):
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
    mesh_json = submit_sparql_query(query_body)
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
    tns1 = get_mesh_tree_numbers(mesh_id1)
    tns2 = get_mesh_tree_numbers(mesh_id2)
    for t1, t2 in itertools.product(tns1, tns2):
        if t1.startswith(t2):
            return True
    return False


def mesh_isa_web(mesh_id1, mesh_id2):
    query_body = """
        SELECT DISTINCT ?o
        FROM <http://id.nlm.nih.gov/mesh>
        WHERE {
          mesh:%s meshv:broaderDescriptor+ ?o .
        }
        """ % mesh_id1
    mesh_json = submit_sparql_query(query_body)
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


def get_mesh_tree_numbers(mesh_id):
    """Return MeSH tree IDs associated with a MeSH ID from the resource file.

    Parameters
    ----------
    mesh_id : str
        The MeSH ID whose tree IDs should be returned.

    Returns
    -------
    list[str]
        A list of MeSH tree IDs.
    """
    return mesh_id_to_tree_numbers.get(mesh_id, [])


def get_mesh_tree_numbers_from_web(mesh_id):
    """Return MeSH tree IDs associated with a MeSH ID from the web.

    Parameters
    ----------
    mesh_id : str
        The MeSH ID whose tree IDs should be returned.

    Returns
    -------
    list[str]
        A list of MeSH tree IDs.
    """
    query_body = """
        SELECT DISTINCT ?tn
        FROM <http://id.nlm.nih.gov/mesh>
        WHERE {
          mesh:%s meshv:treeNumber ?tn
        }
        """ % mesh_id
    mesh_json = submit_sparql_query(query_body)
    if mesh_json is None:
        return []
    try:
        tree_numbers = []
        results = mesh_json['results']['bindings']
        for res in results:
            tree_uri = res['tn']['value']
            m = re.match('http://id.nlm.nih.gov/mesh/([A-Z0-9.]*)', tree_uri)
            tree = m.groups()[0]
            tree_numbers.append(tree)
        return tree_numbers
    except Exception:
        return []


def has_tree_prefix(mesh_id, tree_prefix):
    """Return True if the given MeSH ID has the given tree prefix."""
    tree_numbers = get_mesh_tree_numbers(mesh_id)
    return any(tn.startswith(tree_prefix) for tn in tree_numbers)


def is_disease(mesh_id):
    """Return True if the given MeSH ID is a disease."""
    return has_tree_prefix(mesh_id, 'C')


def is_molecular(mesh_id):
    """Return True if the given MeSH ID is a chemical or drug (incl protein)."""
    return has_tree_prefix(mesh_id, 'D')


def is_enzyme(mesh_id):
    """Return True if the given MeSH ID is an enzyme."""
    return has_tree_prefix(mesh_id, 'D08')


def is_protein(mesh_id):
    """Return True if the given MeSH ID is a protein."""
    return has_tree_prefix(mesh_id, 'D12')


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
    res = get_db_mapping(mesh_id)
    if res and res[0] == 'GO':
        return res[1]
    return None


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
    return get_mesh_id_from_db_id('GO', go_id)


def get_db_mapping(mesh_id):
    """Return mapping to another name space for a MeSH ID, if it exists.

    Parameters
    ----------
    mesh_id : str
        The MeSH ID whose mappings is to be returned.

    Returns
    -------
    tuple or None
        A tuple consisting of a DB namespace and ID for the mapping or None
        if not available.
    """
    return mesh_to_db.get(mesh_id)


def get_mesh_id_from_db_id(db_ns, db_id):
    """Return a MeSH ID mapped from another namespace and ID.

    Parameters
    ----------
    db_ns : str
        A namespace corresponding to db_id.
    db_id : str
        An ID in the given namespace.

    Returns
    -------
    str or None
        The MeSH ID corresponding to the given namespace and ID if available,
        otherwise None.
    """
    return db_to_mesh.get((db_ns, db_id))


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
