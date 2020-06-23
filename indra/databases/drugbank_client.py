import os
from indra.util import read_unicode_csv


mappings_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             os.pardir, 'resources', 'drugbank_mappings.tsv')


def get_chebi_id(drugbank_id):
    return drugbank_chebi.get(drugbank_id)


def get_chembl_id(drugbank_id):
    return drugbank_chembl.get(drugbank_id)


def get_drugbank_from_chebi(chebi_id):
    return chebi_drugbank.get(chebi_id)


def get_drugbank_from_chembl(chembl_id):
    return chembl_drugbank.get(chembl_id)


def load_mappings():
    drugbank_chebi = {}
    chebi_drugbank = {}
    drugbank_chembl = {}
    chembl_drugbank = {}
    for row in read_unicode_csv(mappings_file, delimiter='\t', skiprows=1):
        drugbank_id, db_ns, db_id, source = row
        if db_ns == 'CHEBI':
            chebi_id = 'CHEBI:%s' % db_id
            if drugbank_id in drugbank_chebi or chebi_id in chebi_drugbank:
                import ipdb; ipdb.set_trace()
            drugbank_chebi[drugbank_id] = chebi_id
            chebi_drugbank[chebi_id] = drugbank_id
        elif db_ns == 'CHEMBL':
            if drugbank_id in drugbank_chembl or db_id in chembl_drugbank:
                import ipdb; ipdb.set_trace()
            drugbank_chembl[drugbank_id] = db_id
            chembl_drugbank[db_id] = drugbank_id
    return drugbank_chebi, chebi_drugbank, drugbank_chembl, chembl_drugbank


drugbank_chebi, chebi_drugbank, drugbank_chembl, chembl_drugbank = \
    load_mappings()
