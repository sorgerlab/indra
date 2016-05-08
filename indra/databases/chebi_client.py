import os
import csv
from functools32 import lru_cache

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

def get_pubchem_id(chebi_id):
    return chebi_pubchem.get(chebi_id)
