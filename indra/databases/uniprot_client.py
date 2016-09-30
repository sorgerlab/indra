from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import csv
import rdflib
import logging
import requests
try:
    # Python 3
    from functools import lru_cache
    from urllib.error import HTTPError
except ImportError:
    # Python 2
    from functools32 import lru_cache
    from urllib2 import HTTPError
from indra.util import read_unicode_csv

logger = logging.getLogger('uniprot')

uniprot_url = 'http://www.uniprot.org/uniprot/'

rdf_prefixes = """
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX db: <http://purl.uniprot.org/database/>
    PREFIX faldo: <http://biohackathon.org/resource/faldo#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> """


@lru_cache(maxsize=10000)
def query_protein(protein_id):
    """Return the UniProt entry as an RDF graph for the given UniProt ID.

    Parameters
    ----------
    protein_id : str
        UniProt ID to be queried.

    Returns
    -------
    g : rdflib.Graph
        The RDF graph corresponding to the UniProt entry.
    """
    # Try looking up a primary ID if the given one
    # is a secondary ID
    try:
        prim_ids = uniprot_sec[protein_id]
        protein_id = prim_ids[0]
    except KeyError:
        pass
    url = uniprot_url + protein_id + '.rdf'
    g = rdflib.Graph()
    try:
        g.parse(url)
    except HTTPError:
        logger.warning('Could not find protein with id %s' % protein_id)
        return None
    # Check if the entry has been replaced by a new entry
    query = rdf_prefixes + """
        SELECT ?res2
        WHERE {
            ?res1 up:replacedBy ?res2 .
            }
        """
    res = g.query(query)
    if res:
        term = [r for r in res][0][0]
        replaced_by_id = term.split('/')[-1]
        return query_protein(replaced_by_id)
    return g

def get_family_members(family_name, human_only=True):
    """Return the HGNC gene symbols which are the members of a given family.

    Parameters
    ----------
    family_name : str
        Family name to be queried.
    human_only : bool
        If True, only human proteins in the family will be returned.
        Default: True

    Returns
    -------
    gene_names : list
        The HGNC gene symbols corresponding to the given family.
    """
    data = {'query': 'family:%s' % family_name, 
            'format': 'list'}
    if human_only:
        data['fil'] = 'organism:human'
    res = requests.get(uniprot_url, params=data)
    if not res.status_code == 200 or not res.text:
        return None
    # res.text gets us the Unicode
    html = res.text
    protein_list = html.strip().split('\n')
    gene_names = []
    for p in protein_list:
        gene_name = get_gene_name(p)
        gene_names.append(gene_name)
    return gene_names

def get_mnemonic(protein_id, web_fallback=True):
    """Return the UniProt mnemonic for the given UniProt ID.

    Parameters
    ----------
    protein_id : str
        UniProt ID to be mapped.
    web_fallback : Optional[bool]
        If True and the offline lookup fails, the UniProt web service
        is used to do the query.

    Returns
    -------
    mnemonic : str
        The UniProt mnemonic corresponding to the given Uniprot ID.
    """
    try:
        mnemonic = uniprot_mnemonic[protein_id]
        return mnemonic
    except KeyError:
        pass
    if not web_fallback:
        return None
    g = query_protein(protein_id)
    if g is None:
        return None
    query = rdf_prefixes + """
        SELECT ?mnemonic
        WHERE {
            ?r up:mnemonic ?mnemonic .
        }
        """
    res = g.query(query)
    if res:
        mnemonic = [r for r in res][0][0].toPython()
        return mnemonic
    else:
        return None

def get_id_from_mnemonic(uniprot_mnemonic):
    """Return the UniProt ID for the given UniProt mnemonic.

    Parameters
    ----------
    uniprot_mnemonic : str
        UniProt mnemonic to be mapped.

    Returns
    -------
    uniprot_id : str
        The UniProt ID corresponding to the given Uniprot mnemonic.
    """
    try:
        uniprot_id = uniprot_mnemonic_reverse[uniprot_mnemonic]
        return uniprot_id
    except KeyError:
        return None

def get_gene_name(protein_id, web_fallback=True):
    """Return the gene name for the given UniProt ID.

    This is an alternative to get_hgnc_name and is useful when
    HGNC name is not availabe (for instance, when the organism
    is not homo sapiens).

    Parameters
    ----------
    protein_id : str
        UniProt ID to be mapped.
    web_fallback : Optional[bool]
        If True and the offline lookup fails, the UniProt web service
        is used to do the query.

    Returns
    -------
    gene_name : str
        The gene name corresponding to the given Uniprot ID.
    """
    try:
        gene_name = uniprot_gene_name[protein_id]
        # Handle empty string
        if not gene_name:
            return None
        return gene_name
    except KeyError:
        pass
    if not web_fallback:
        return None

    g = query_protein(protein_id)
    if g is None:
        return None
    query = rdf_prefixes + """
        SELECT ?name
        WHERE {
            ?gene a up:Gene .
            ?gene skos:prefLabel ?name .
            }
        """
    res = g.query(query)
    if res:
        gene_name = [r for r in res][0][0].toPython()
        if not gene_name:
            return None
        return gene_name
    return None

@lru_cache(maxsize=1000)
def get_sequence(protein_id):
    try:
        prim_ids = uniprot_sec[protein_id]
        protein_id = prim_ids[0]
    except KeyError:
        pass
    url = uniprot_url + '%s.fasta' % protein_id
    res = requests.get(url)
    if not res.status_code == 200:
        logger.warning('Could not find sequence for protein %s' % protein_id)
        return None
    # res.text is Unicode
    lines = res.text.splitlines()
    seq = (''.join(lines[1:])).replace('\n','')
    return seq


def get_modifications(protein_id):
    g = query_protein(protein_id)
    if g is None:
        return None
    query = rdf_prefixes + """
        SELECT ?beg_pos ?comment
        WHERE {
            ?mod_res a up:Modified_Residue_Annotation .
            ?mod_res rdfs:comment ?comment .
            ?mod_res up:range ?range .
            ?range faldo:begin ?beg .
            ?range faldo:end ?end .
            ?beg a faldo:ExactPosition .
            ?beg faldo:position ?beg_pos .
            FILTER (?beg = ?end)
            }
        """
    res = g.query(query)
    mods = []
    for r in res:
        mod_pos = r[0].value
        # "Phosphothreonine; by autocatalysis"
        # "Phosphothreonine; by MAP2K1 and MAP2K2"
        # TODO: take into account the comment after the ;?
        mod_res = r[1].value.split(';')[0]
        mods.append((mod_res, mod_pos))
    return mods


def verify_location(protein_id, residue, location):
    """Return True if the residue is at the given location in the UP sequence.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the protein whose sequence is used as reference.
    residue : str
        A single character amino acid symbol (Y, S, T, V, etc.)
    location : str
        The location on the protein sequence (starting at 1) at which the
        residue should be checked against the reference sequence.

    Returns
    -------
    True if the given residue is at the given position in the sequence
    corresponding to the given UniProt ID, otherwise False.
    """
    seq = get_sequence(protein_id)
    # If we couldn't get the sequence (can happen due to web service hiccups)
    # don't throw the statement away by default
    if seq is None:
        return True
    try:
        loc_int = int(location)
    except ValueError:
        logger.warning('Invalid location %s' % location)
        loc_int = -1

    if (loc_int < 1) or (loc_int > len(seq)):
        return False
    elif seq[loc_int - 1] == residue:
        return True
    return False


def verify_modification(protein_id, residue, location=None):
    """Return True if the residue at the given location has a known modifiation. 

    Parameters
    ----------
    protein_id : str
        UniProt ID of the protein whose sequence is used as reference.
    residue : str
        A single character amino acid symbol (Y, S, T, V, etc.)
    location : Optional[str]
        The location on the protein sequence (starting at 1) at which the
        modification is checked.

    Returns
    -------
    True if the given residue is reported to be modified at the given position 
    in the sequence corresponding to the given UniProt ID, otherwise False.
    If location is not given, we only check if there is any residue of the
    given type that is modified.
    """
    mods = get_modifications(protein_id)
    mod_locs = [m[1] for m in mods]
    seq = get_sequence(protein_id)
    if location:
        if not verify_location(protein_id, residue, location):
            return False
        try:
            mod_idx = mod_locs.index(location)
        except ValueError:
            return False
        return True
    else:
        for ml in mod_locs:
            if seq[ml - 1] == residue:
                return True
        return False

def _build_uniprot_entries():
    up_entries_file = os.path.dirname(os.path.abspath(__file__)) + \
        '/../resources/uniprot_entries.tsv'
    uniprot_gene_name = {}
    uniprot_mnemonic = {}
    uniprot_mnemonic_reverse = {}
    try:
        csv_rows = read_unicode_csv(up_entries_file, delimiter='\t')
        # Skip the header row
        next(csv_rows)
        for row in csv_rows:
            up_id = row[0]
            gene_name = row[1]
            up_mnemonic = row[3]
            uniprot_gene_name[up_id] = gene_name
            uniprot_mnemonic[up_id] = up_mnemonic
            uniprot_mnemonic_reverse[up_mnemonic] = up_id
    except IOError:
        pass
    return uniprot_gene_name, uniprot_mnemonic, uniprot_mnemonic_reverse

def _build_uniprot_hgnc():
    hgnc_file = os.path.dirname(os.path.abspath(__file__)) +\
                '/../resources/hgnc_entries.txt'
    try:
        csv_rows = read_unicode_csv(hgnc_file, delimiter='\t')
        # Skip the header row
        next(csv_rows)
        uniprot_hgnc = {}
        for row in csv_rows:
            hgnc_name = row[1]
            uniprot_id = row[6]
            if uniprot_id:
                uniprot_hgnc[uniprot_id] = hgnc_name
    except IOError:
        uniprot_hgnc = {}
    return uniprot_hgnc

def _build_uniprot_sec():
    # File containing secondary accession numbers mapped
    # to primary accession numbers
    sec_file = os.path.dirname(os.path.abspath(__file__)) +\
                '/../resources/uniprot_sec_ac.txt'
    try:
        uniprot_sec = {}
        lines = open(sec_file, 'rt').readlines()
        for i, l in enumerate(lines):
            if l.startswith('Secondary AC'):
                entry_lines = lines[i+2:]

        for l in entry_lines:
            sec_id, prim_id = l.split()
            try:
                uniprot_sec[sec_id].append(prim_id)
            except KeyError:
                uniprot_sec[sec_id] = [prim_id]
    except IOError:
        uniprot_sec = {}
    return uniprot_sec

def _build_uniprot_subcell_loc():
    fname = os.path.dirname(os.path.abspath(__file__)) +\
                '/../resources/uniprot_subcell_loc.tsv'
    try:
        csv_rows = read_unicode_csv(fname, delimiter='\t')
        # Skip the header row
        next(csv_rows)
        subcell_loc = {}
        for row in csv_rows:
            loc_id = row[0]
            loc_alias = row[3]
            subcell_loc[loc_id] = loc_alias
    except IOError:
        subcell_loc = {}
    return subcell_loc

uniprot_gene_name, uniprot_mnemonic, uniprot_mnemonic_reverse = \
    _build_uniprot_entries()
uniprot_sec = _build_uniprot_sec()
uniprot_subcell_loc = _build_uniprot_subcell_loc()
