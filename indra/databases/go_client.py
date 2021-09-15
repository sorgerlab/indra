"""A client to the Gene Ontology."""
import re
import logging
from typing import Union
from indra.databases.obo_client import OboClient

logger = logging.getLogger(__name__)

_client = OboClient(prefix='go')


def get_go_label(go_id):
    """Get label corresponding to a given GO identifier.

    Parameters
    ----------
    go_id : str
        The GO identifier. Should include the `GO:` prefix, e.g., `GO:1903793`
        (positive regulation of anion transport).

    Returns
    -------
    str
        Label corresponding to the GO ID.
    """
    return _client.get_name_from_id(go_id)


def get_go_id_from_label(label):
    """Get ID corresponding to a given GO label.

    Parameters
    ----------
    label : str
        The GO label to get the ID for.

    Returns
    -------
    str
        Identifier corresponding to the GO label, starts with GO:.
    """
    return _client.get_id_from_name(label)


def get_go_id_from_label_or_synonym(label):
    """Get ID corresponding to a given GO label or synonym

    Parameters
    ----------
    label : str
        The GO label or synonym to get the ID for.

    Returns
    -------
    str
        Identifier corresponding to the GO label or synonym, starts with GO:.
    """
    return _client.get_id_from_name_or_synonym(label)


def get_primary_id(go_id):
    """Get primary ID corresponding to an alternative/deprecated GO ID.

    Parameters
    ----------
    go_id : str
        The GO ID to get the primary ID for.

    Returns
    -------
    str
        Primary identifier corresponding to the given ID.
    """
    return _client.get_id_from_alt_id(go_id)


def get_valid_location(loc):
    """Return a valid GO label based on an ID, label or synonym.

    The rationale behind this function is that many sources produce
    cellular locations that are arbitrarily either GO IDs (sometimes
    without the prefix and sometimes outdated) or labels or synonyms.
    This function handles all these cases and returns a valid GO label
    in case one is available, otherwise None.

    Parameters
    ----------
    loc : txt
        The location that needst o be canonicalized.

    Returns
    -------
    str or None
        The valid location string is available, otherwise None.
    """
    if not loc:
        return None
    # If it's actually a GO ID, we do some validation and use it. If it is
    # a text label then we look up the GO ID for it
    if re.match(r'^(GO:)?\d+$', loc):
        if not loc.startswith('GO:'):
            loc = 'GO:' + loc
        go_id = loc
        prim_id = get_primary_id(go_id)
        if prim_id:
            go_id = prim_id
    else:
        go_id = get_go_id_from_label_or_synonym(loc)
        if not go_id:
            return None
    # If we managed to get a GO ID either way, we get its label and return it
    # with some extra caution to not return a None name under any
    # circumstances
    if go_id:
        loc = get_go_label(go_id)
        if loc:
            return loc
    return None


def get_namespace(go_id: str) -> Union[str, None]:
    """Return the GO namespace associated with a GO ID.

    Parameters
    ----------
    go_id :
        The GO ID to get the namespace for

    Returns
    -------
    :
        The GO namespace for the given ID. This is one of
        molecular_function, biological_process or cellular_component.
        If the GO ID is not available as an entry, None is returned.
    """
    return _client.entries.get(go_id, {}).get('namespace')