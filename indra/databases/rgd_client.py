"""A client for accessing RGD rat gene data."""

from collections import defaultdict
from typing import List, Union

from indra.util import read_unicode_csv
from indra.resources import get_resource_path


def get_id_from_name(name: str) -> Union[str, None]:
    """Return an RGD ID from an RGD gene symbol.

    Parameters
    ----------
    name :
        The RGD gene symbol whose ID will be returned.

    Returns
    -------
    :
        The RGD ID (without prefix) or None if not available.
    """
    return rgd_name_to_id.get(name)


def get_name_from_id(rgd_id: str) -> Union[str, None]:
    """Return the RGD gene symbol for a given RGD ID.

    Parameters
    ----------
    rgd_id :
        The RGD ID (without prefix) whose symbol will be returned.

    Returns
    -------
    :
        The RGD symbol for the given ID or None if not available.
    """
    if rgd_id and rgd_id.startswith('RGD:'):
        rgd_id = rgd_id[4:]
    return rgd_id_to_name.get(rgd_id)


def get_synonyms(rgd_id: str) -> List[str]:
    """Return the synonyms for an RGD ID.

    Parameters
    ----------
    rgd_id :
        An RGD ID, without prefix.

    Returns
    -------
    :
        The list of synonyms corresponding to the RGD ID, or an empty list
        if not available.
    """
    if rgd_id and rgd_id.startswith('RGD:'):
        rgd_id = rgd_id[4:]
    return rgd_synonyms.get(rgd_id, [])


def get_id_from_name_synonym(name_synonym: str) -> Union[None, str, List[str]]:
    """Return an RGD ID from an RGD gene symbol or synonym.

    If the given name or synonym is the official symbol of a gene, its
    ID is returned. If the input is a synonym, it can correspond to
    one or more genes. If there is a single gene whose synonym matches
    the input, the ID is returned as a string. If multiple genes share
    the given synonym, their IDs are returned in a list. If the input
    doesn't match any names or synonyms, None is returned.

    Parameters
    ----------
    name_synonym :
        The RGD gene symbol or synonym whose ID will be returned.

    Returns
    -------
    :
        The RGD ID (without prefix) of a single gene, a list of RGD IDs,
        or None.
    """
    rgd_id = rgd_name_to_id.get(name_synonym)
    if rgd_id:
        return rgd_id
    rgd_ids = rgd_synonyms_reverse.get(name_synonym)
    if rgd_ids:
        if len(rgd_ids) == 1:
            return rgd_ids[0]
        else:
            return rgd_ids
    return None


def get_ensembl_id(rgd_id: str) -> Union[str, None]:
    """Return the Ensembl ID for an RGD ID.

    Parameters
    ----------
    rgd_id :
        An RGD ID, without prefix.

    Returns
    -------
    :
        A list of Ensembl IDs corresponding to the RGD ID,
        or None if not available.
    """
    return rgd_id_to_ensembl.get(rgd_id)


def _read_rgd():
    fname = get_resource_path('rgd_entries.tsv')
    rgd_id_to_name = {}
    rgd_name_to_id = {}
    rgd_synonyms = {}
    rgd_synonyms_reverse = defaultdict(list)
    rgd_id_to_ensembl = {}
    ensemble_id_to_rgd = {}
    for rgd_id, name, synonyms_str, ensembl_id in \
            read_unicode_csv(fname, '\t'):
        if name:
            rgd_id_to_name[rgd_id] = name
            rgd_name_to_id[name] = rgd_id
        if synonyms_str:
            synonyms = synonyms_str.split(';')
            rgd_synonyms[rgd_id] = synonyms
            for synonym in synonyms:
                rgd_synonyms_reverse[synonym].append(rgd_id)
        if ensembl_id:
            ensemble_ids = ensembl_id.split(';')
            rgd_id_to_ensembl[rgd_id] = ensemble_ids
            for ensemble_id in ensemble_ids:
                ensemble_id_to_rgd[ensemble_id] = rgd_id

    return rgd_id_to_name, rgd_name_to_id, rgd_synonyms, \
        dict(rgd_synonyms_reverse), rgd_id_to_ensembl, ensemble_id_to_rgd


rgd_id_to_name, rgd_name_to_id, rgd_synonyms, rgd_synonyms_reverse, \
   rgd_id_to_ensembl, ensemble_id_to_rgd = _read_rgd()
