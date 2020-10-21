"""Client for interacting with DrugBank entries."""
import os
from indra.util import read_unicode_csv


def get_db_mapping(drugbank_id, db_ns):
    """Return a mapping for a DrugBank ID to a given name space.

    Parameters
    ----------
    drugbank_id : str
        DrugBank ID to map.
    db_ns : str
        The database name space to map to.

    Returns
    -------
    str or None
        The ID mapped to the given name space or None if not available.
    """
    return drugbank_to_db.get((drugbank_id, db_ns))


def get_drugbank_id_from_db_id(db_ns, db_id):
    """Return DrugBank ID from a database name space and ID.

    Parameters
    ----------
    db_ns : str
        Database name space.
    db_id : str
        Database ID.

    Returns
    -------
    str or None
        The mapped DrugBank ID or None if not available.
    """
    return db_to_drugbank.get((db_ns, db_id))


def get_chebi_id(drugbank_id):
    """Return a mapping for a DrugBank ID to CHEBI.

    Parameters
    ----------
    drugbank_id : str
        DrugBank ID to map.

    Returns
    -------
    str or None
        The ID mapped to CHEBI or None if not available.
    """
    return get_db_mapping(drugbank_id, 'CHEBI')


def get_chembl_id(drugbank_id):
    """Return a mapping for a DrugBank ID to CHEMBL.

    Parameters
    ----------
    drugbank_id : str
        DrugBank ID to map.

    Returns
    -------
    str or None
        The ID mapped to CHEMBL or None if not available.
    """
    return get_db_mapping(drugbank_id, 'CHEMBL')


def get_drugbank_id_from_chebi_id(chebi_id):
    """Return DrugBank ID from a CHEBI ID.

    Parameters
    ----------
    chebi_id : str
        CHEBI ID to map.

    Returns
    -------
    str or None
        The mapped DrugBank ID or None if not available.
    """
    return get_drugbank_id_from_db_id('CHEBI', chebi_id)


def get_drugbank_id_from_chembl_id(chembl_id):
    """Return DrugBank ID from a CHEMBL ID.

    Parameters
    ----------
    chembl_id : str
        CHEMBL ID to map.

    Returns
    -------
    str or None
        The mapped DrugBank ID or None if not available.
    """
    return get_drugbank_id_from_db_id('CHEMBL', chembl_id)


def get_drugbank_name(drugbank_id):
    """Return the DrugBank standard name for a given DrugBank ID.

    Parameters
    ----------
    drugbank_id : str
        DrugBank ID to get the name for

    Returns
    -------
    str or None
        The name corresponding to the given DrugBank ID or None if not
        available.
    """
    return drugbank_names.get(drugbank_id)


mappings_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             os.pardir, 'resources', 'drugbank_mappings.tsv')


def _load_mappings():
    drugbank_to_db = {}
    db_to_drugbank = {}
    drugbank_names = {}
    to_db_ambigs = set()
    db_to_ambigs = set()
    for drugbank_id, db_ns, db_id, source in \
            read_unicode_csv(mappings_file, delimiter='\t', skiprows=1):
        # We skip DBSALTs for now, see https://github.com/pyobo/pyobo/issues/80
        if drugbank_id.startswith('DBSALT'):
            continue
        if db_ns == 'CHEBI':
            db_id = 'CHEBI:%s' % db_id
        if db_ns == 'NAME':
            drugbank_names[drugbank_id] = db_id
            continue
        key = (db_ns, db_id)
        if key in db_to_drugbank and db_to_drugbank[key] != drugbank_id:
            db_to_ambigs.add(key)
        db_to_drugbank[key] = drugbank_id
        key = (drugbank_id, db_ns)
        if key in drugbank_to_db and \
                drugbank_to_db[key] != db_id:
            to_db_ambigs.add(key)
        drugbank_to_db[key] = db_id
    db_to_drugbank = {k: v for k, v in db_to_drugbank.items()
                      if k not in db_to_ambigs}
    drugbank_to_db = {k: v for k, v in drugbank_to_db.items()
                      if k not in to_db_ambigs}
    return drugbank_to_db, db_to_drugbank, drugbank_names


drugbank_to_db, db_to_drugbank, drugbank_names = _load_mappings()
