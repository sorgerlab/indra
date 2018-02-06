from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import dirname, abspath, join
try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache
from indra.util import read_unicode_csv


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
    """Return a ChEMBL ID from a given ChEBI ID.

    Parameters
    ----------
    chebi_id : str
        ChEBI ID to be converted.

    Returns
    -------
    chembl_id : str
        ChEMBL ID corresponding to the given ChEBI ID. If the lookup fails,
        None is returned.
    """
    return chebi_chembl.get(chebi_id)


def get_chebi_id_from_cas(cas_id):
    """Return a ChEBI ID corresponding to the given CAS ID.

    Parameters
    ----------
    cas_id : str
        The CAS ID to be converted.

    Parameters
    ----------
    chebi_id : str
        The ChEBI ID corresponding to the given CAS ID. If the lookup
        fails, None is returned.
    """
    return cas_chebi.get(cas_id)


def _read_chebi_to_pubchem():
    chebi_to_pubchem_file = join(dirname(abspath(__file__)),
                                 '../resources/chebi_to_pubchem.tsv')
    csv_reader = read_unicode_csv(chebi_to_pubchem_file, delimiter='\t')
    chebi_pubchem = {}
    pubchem_chebi = {}
    for row in csv_reader:
        chebi_pubchem[row[0]] = row[1]
        pubchem_chebi[row[1]] = row[0]
    return chebi_pubchem, pubchem_chebi


def _read_chebi_to_chembl():
    chebi_to_chembl_file = join(dirname(abspath(__file__)),
                                '../resources/chebi_to_chembl.tsv')
    csv_reader = read_unicode_csv(chebi_to_chembl_file, delimiter='\t')
    chebi_chembl = {}
    for row in csv_reader:
        chebi_chembl[row[0]] = row[1]
    return chebi_chembl


def _read_cas_to_chebi():
    cas_to_chebi_file = join(dirname(abspath(__file__)),
                             '../resources/cas_to_chebi.tsv')
    csv_reader = read_unicode_csv(cas_to_chebi_file, delimiter='\t')
    cas_chebi = {}
    next(csv_reader)
    for row in csv_reader:
        cas_chebi[row[0]] = row[1]
    # These are missing from the resource but appear often, so we map
    # them manually
    extra_entries = {'24696-26-2': '17761',
                     '23261-20-3': '18035',
                     '165689-82-7': '16618'}
    cas_chebi.update(extra_entries)
    return cas_chebi


chebi_pubchem, pubchem_chebi = _read_chebi_to_pubchem()
chebi_chembl = _read_chebi_to_chembl()
cas_chebi = _read_cas_to_chebi()
