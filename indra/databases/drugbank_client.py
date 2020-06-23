import os
from indra.util import read_unicode_csv


mappings_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             os.pardir, 'resources', 'drugbank_mappings.tsv')


def get_db_mapping(drugbank_id):
    return drugbank_to_db.get(drugbank_id)


def get_drugbank_id_from_db_id(db_ns, db_id):
    return db_to_drugbank.get((db_ns, db_id))


def get_chebi_id(drugbank_id):
    res = get_db_mapping(drugbank_id)
    if res and res[0] == 'CHEBI':
        return res[1]
    return None


def get_chembl_id(drugbank_id):
    res = get_db_mapping(drugbank_id)
    if res and res[0] == 'CHEMBL':
        return res[1]
    return None


def get_drugbank_id_from_chebi_id(chebi_id):
    return get_drugbank_id_from_db_id('CHEBI', chebi_id)


def get_drugbank_id_from_chembl_id(chembl_id):
    return get_drugbank_id_from_db_id('CHEMBL', chembl_id)


def load_mappings():
    drugbank_to_db = {}
    db_to_drugbank = {}
    to_db_ambigs = set()
    db_to_ambigs = set()
    for drugbank_id, db_ns, db_id, source in \
            read_unicode_csv(mappings_file, delimiter='\t', skiprows=1):
        if db_ns == 'CHEBI':
            db_id = 'CHEBI:%s' % db_id
        db_key = (db_ns, db_id)
        if db_key in db_to_drugbank and db_to_drugbank[db_key] != drugbank_id:
            db_to_ambigs.add(db_key)
        if drugbank_id in drugbank_to_db and \
                drugbank_to_db[drugbank_id] != db_key:
            to_db_ambigs.add(drugbank_id)
        db_to_drugbank[db_key]  = drugbank_id
        drugbank_to_db[drugbank_id] = db_key
    db_to_drugbank = {k: v for k, v in db_to_drugbank.items()
                      if k not in db_to_ambigs}
    drugbank_to_db = {k: v for k, v in drugbank_to_db.items()
                      if k not in to_db_ambigs}
    return drugbank_to_db, db_to_drugbank


drugbank_to_db, db_to_drugbank = load_mappings()
