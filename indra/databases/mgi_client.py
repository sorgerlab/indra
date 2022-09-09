"""A client for accessing MGI mouse gene data."""

from collections import defaultdict
from typing import List, Union

from indra.util import read_unicode_csv
from indra.resources import get_resource_path


def get_id_from_name(name: str) -> Union[str, None]:
    """Return an MGI ID from an MGI gene symbol.

    Parameters
    ----------
    name :
        The MGI gene symbol whose ID will be returned.

    Returns
    -------
    :
        The MGI ID (without prefix) or None if not available.
    """
    return mgi_name_to_id.get(name)


def get_name_from_id(mgi_id: str) -> Union[str, None]:
    """Return the MGI gene symbol for a given MGI ID.

    Parameters
    ----------
    mgi_id :
        The MGI ID (without prefix) whose symbol will be returned.

    Returns
    -------
    :
        The MGI symbol for the given ID or None if not available.
    """
    if mgi_id and mgi_id.startswith('MGI:'):
        mgi_id = mgi_id[4:]
    return mgi_id_to_name.get(mgi_id)


def get_synonyms(mgi_id: str) -> List[str]:
    """Return the synonyms for an MGI ID.

    Parameters
    ----------
    mgi_id :
        An MGI ID, without prefix.

    Returns
    -------
    :
        The list of synonyms corresponding to the MGI ID, or an empty list
        if not available.
    """
    if mgi_id and mgi_id.startswith('MGI:'):
        mgi_id = mgi_id[4:]
    return mgi_synonyms.get(mgi_id, [])


def get_id_from_name_synonym(name_synonym: str) -> Union[None, str, List[str]]:
    """Return an MGI ID from an MGI gene symbol or synonym.

    If the given name or synonym is the official symbol of a gene, its
    ID is returned. If the input is a synonym, it can correspond to
    one or more genes. If there is a single gene whose synonym matches
    the input, the ID is returned as a string. If multiple genes share
    the given synonym, their IDs are returned in a list. If the input
    doesn't match any names or synonyms, None is returned.

    Parameters
    ----------
    name_synonym :
        The MGI gene symbol or synonym whose ID will be returned.

    Returns
    -------
    :
        The MGI ID (without prefix) of a single gene, a list of MGI IDs,
        or None.
    """
    mgi_id = mgi_name_to_id.get(name_synonym)
    if mgi_id:
        return mgi_id
    mgi_ids = mgi_synonyms_reverse.get(name_synonym)
    if mgi_ids:
        if len(mgi_ids) == 1:
            return mgi_ids[0]
        else:
            return mgi_ids
    return None


def _read_mgi():
    fname = get_resource_path('mgi_entries.tsv')
    mgi_id_to_name = {}
    mgi_name_to_id = {}
    mgi_synonyms = {}
    mgi_synonyms_reverse = defaultdict(list)
    for mgi_id, name, synonyms_str in read_unicode_csv(fname, '\t'):
        if name:
            mgi_id_to_name[mgi_id] = name
            mgi_name_to_id[name] = mgi_id
        if synonyms_str:
            synonyms = synonyms_str.split('|')
            mgi_synonyms[mgi_id] = synonyms
            for synonym in synonyms:
                mgi_synonyms_reverse[synonym].append(mgi_id)

    return mgi_id_to_name, mgi_name_to_id, mgi_synonyms, \
        dict(mgi_synonyms_reverse)


mgi_id_to_name, mgi_name_to_id, mgi_synonyms, mgi_synonyms_reverse = _read_mgi()
