from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import dirname, abspath, join
try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache
from indra.util import read_unicode_csv

def read_chebi_to_pubchem():
    chebi_to_pubchem_file = join(dirname(abspath(__file__)),
                                 '../resources/chebi_to_pubchem.tsv')
    csv_reader = read_unicode_csv(chebi_to_pubchem_file, delimiter='\t')
    chebi_pubchem = {}
    pubchem_chebi = {}
    for row in csv_reader:
        chebi_pubchem[row[0]] = row[1]
        pubchem_chebi[row[1]] = row[0]
    return chebi_pubchem, pubchem_chebi

def read_chebi_to_chembl():
    chebi_to_chembl_file = join(dirname(abspath(__file__)),
                                '../resources/chebi_to_chembl.tsv')
    csv_reader = read_unicode_csv(chebi_to_chembl_file, delimiter='\t')
    chebi_chembl = {}
    for row in csv_reader:
        chebi_chembl[row[0]] = row[1]
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
