import os
import csv
from functools32 import lru_cache

def read_chebi_to_pubchem():
    # Based on ftp://ftp.ebi.ac.uk/pub/databases/chebi/
    #                Flat_file_tab_delimited/reference.tsv.gz
    chebi_to_pubchem_file = os.path.dirname(os.path.abspath(__file__)) + \
                            '/../resources/chebi_to_pubchem.tsv'
    try:
        fh = open(chebi_to_pubchem_file, 'rt')
        rd = csv.reader(fh, delimiter='\t')
        chebi_pubchem = {}
        for row in rd:
            chebi_pubchem[row[0]] = row[1]
    except IOError:
        chebi_pubchem = {}
    return chebi_pubchem

def read_chebi_to_chembl():
    # Based on ftp://ftp.ebi.ac.uk/pub/databases/chebi/
    #                Flat_file_tab_delimited/reference.tsv.gz
    chebi_to_chembl_file = os.path.dirname(os.path.abspath(__file__)) + \
                            '/../resources/chebi_to_chembl.tsv'
    try:
        fh = open(chebi_to_chembl_file, 'rt')
        rd = csv.reader(fh, delimiter='\t')
        chebi_chembl = {}
        for row in rd:
            chebi_chembl[row[0]] = row[1]
    except IOError:
        chebi_chembl = {}
    return chebi_chembl

chebi_pubchem = read_chebi_to_pubchem()
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
    return chebi_pubchem.get(chebi_id)

def get_chembl_id(chebi_id):
    return chebi_chembl.get(chebi_id)
