from os.path import dirname, abspath, join

import requests
import logging

from functools import lru_cache

from lxml import etree, objectify

from indra.util import read_unicode_csv


logger = logging.getLogger(__name__)


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

    Returns
    -------
    chebi_id : str
        The ChEBI ID corresponding to the given CAS ID. If the lookup
        fails, None is returned.
    """
    return cas_chebi.get(cas_id)


def get_chebi_name_from_id(chebi_id, offline=False):
    """Return a ChEBI name corresponding to the given ChEBI ID.

    Parameters
    ----------
    chebi_id : str
        The ChEBI ID whose name is to be returned.
    offline : Optional[bool]
        Choose whether to allow an online lookup if the local lookup fails. If
        True, the online lookup is not attempted. Default: False.

    Returns
    -------
    chebi_name : str
        The name corresponding to the given ChEBI ID. If the lookup
        fails, None is returned.
    """
    chebi_name = chebi_id_to_name.get(chebi_id)
    if chebi_name is None and not offline:
        chebi_name = get_chebi_name_from_id_web(chebi_id)
    return chebi_name


def get_chebi_id_from_name(chebi_name):
    """Return a ChEBI ID corresponding to the given ChEBI name.

    Parameters
    ----------
    chebi_name : str
        The ChEBI name whose ID is to be returned.

    Returns
    -------
    chebi_id : str
        The ID corresponding to the given ChEBI name. If the lookup
        fails, None is returned.
    """
    chebi_id = chebi_name_to_id.get(chebi_name)
    return chebi_id


def _read_chebi_to_pubchem():
    csv_reader = _read_relative_csv('../resources/chebi_to_pubchem.tsv')
    chebi_pubchem = {}
    pubchem_chebi = {}
    for row in csv_reader:
        chebi_pubchem[row[0]] = row[1]
        pubchem_chebi[row[1]] = row[0]
    return chebi_pubchem, pubchem_chebi


def _read_chebi_to_chembl():
    csv_reader = _read_relative_csv('../resources/chebi_to_chembl.tsv')
    chebi_chembl = {}
    for row in csv_reader:
        chebi_chembl[row[0]] = row[1]
    return chebi_chembl


def _read_cas_to_chebi():
    csv_reader = _read_relative_csv('../resources/cas_to_chebi.tsv')
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


def _read_chebi_names():
    csv_reader = _read_relative_csv('../resources/chebi_names.tsv')
    next(csv_reader)
    chebi_id_to_name = {}
    chebi_name_to_id = {}
    for row in csv_reader:
        chebi_id, name = row
        chebi_id_to_name[chebi_id] = name
        chebi_name_to_id[name] = chebi_id
    return chebi_id_to_name, chebi_name_to_id


def _read_relative_csv(rel_path):
    file_path = join(dirname(abspath(__file__)), rel_path)
    csv_reader = read_unicode_csv(file_path, delimiter='\t')
    return csv_reader


@lru_cache(maxsize=5000)
def get_chebi_name_from_id_web(chebi_id):
    """Return a ChEBI mame corresponding to a given ChEBI ID using a REST API.

    Parameters
    ----------
    chebi_id : str
        The ChEBI ID whose name is to be returned.

    Returns
    -------
    chebi_name : str
        The name corresponding to the given ChEBI ID. If the lookup
        fails, None is returned.
    """
    url_base = 'http://www.ebi.ac.uk/webservices/chebi/2.0/test/'
    url_fmt = url_base + 'getCompleteEntity?chebiId=%s'
    resp = requests.get(url_fmt % chebi_id)
    if resp.status_code != 200:
        logger.warning("Got bad code form CHEBI client: %s" % resp.status_code)
        return None
    tree = etree.fromstring(resp.content)

    # Get rid of the namespaces.
    # Credit: https://stackoverflow.com/questions/18159221/remove-namespace-and-prefix-from-xml-in-python-using-lxml
    for elem in tree.getiterator():
        if not hasattr(elem.tag, 'find'):
            continue  # (1)
        i = elem.tag.find('}')
        if i >= 0:
            elem.tag = elem.tag[i+1:]
    objectify.deannotate(tree, cleanup_namespaces=True)

    elem = tree.find('Body/getCompleteEntityResponse/return/chebiAsciiName')
    if elem is not None:
        return elem.text
    return None


chebi_pubchem, pubchem_chebi = _read_chebi_to_pubchem()
chebi_chembl = _read_chebi_to_chembl()
cas_chebi = _read_cas_to_chebi()
chebi_id_to_name, chebi_name_to_id = _read_chebi_names()
