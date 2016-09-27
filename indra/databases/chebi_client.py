from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import csv
try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache

def read_chebi_to_pubchem():
    # Based on ftp://ftp.ebi.ac.uk/pub/databases/chebi/
    #                Flat_file_tab_delimited/reference.tsv.gz
    chebi_to_pubchem_file = os.path.dirname(os.path.abspath(__file__)) + \
                            '/../resources/chebi_to_pubchem.tsv'
    try:
        fh = open(chebi_to_pubchem_file, 'rt')
        # In Python 2, the unicode_literal delimiter '\t' will give TypeError
        try:
            rd = csv.reader(fh, delimiter='\t')
        except TypeError:
            rd = csv.reader(fh, delimiter='\t'.encode('utf-8'))
        chebi_pubchem = {}
        pubchem_chebi = {}
        for row in rd:
            chebi_pubchem[row[0]] = row[1]
            pubchem_chebi[row[1]] = row[0]
    except IOError:
        chebi_pubchem = {}
        pubchem_chebi = {}
    return chebi_pubchem, pubchem_chebi

def read_chebi_to_chembl():
    # Based on ftp://ftp.ebi.ac.uk/pub/databases/chebi/
    #                Flat_file_tab_delimited/reference.tsv.gz
    chebi_to_chembl_file = os.path.dirname(os.path.abspath(__file__)) + \
                            '/../resources/chebi_to_chembl.tsv'
    try:
        fh = open(chebi_to_chembl_file, 'rt')
        # In Python 2, the unicode_literal delimiter '\t' will give TypeError
        try:
            rd = csv.reader(fh, delimiter='\t')
        except TypeError:
            rd = csv.reader(fh, delimiter='\t'.encode('utf-8'))
        chebi_chembl = {}
        for row in rd:
            chebi_chembl[row[0]] = row[1]
    except IOError:
        chebi_chembl = {}
    return chebi_chembl

chebi_pubchem, pubchem_chebi = read_chebi_to_pubchem()
chebi_chembl = read_chebi_to_chembl()

def get_pubchem_id(chebi_id):
    """Return the PubChem ID corresponding to a given ChEBI ID.

    Parameters
    ----------
    chebi_id : str
        ChEBI ID to be converted.

    Returns
    -------
    pubchem_id : str
        PubChem ID corresponding to the given ChEBI ID. If the lookup fails,
        None is returned.
    """
    pubchem_id = chebi_pubchem.get(chebi_id)
    return pubchem_id

def get_chebi_id_from_pubchem(pubchem_id):
    """Return the ChEBI ID corresponding to a given Pubchem ID.

    Parameters
    ----------
    pubchem_id : str
        Pubchem ID to be converted.

    Returns
    -------
    chebi_id : str
        ChEBI ID corresponding to the given Pubchem ID. If the lookup fails,
        None is returned.
    """
    chebi_id = pubchem_chebi.get(pubchem_id)
    return chebi_id

def get_chembl_id(chebi_id):
    return chebi_chembl.get(chebi_id)
