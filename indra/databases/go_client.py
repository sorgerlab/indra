"""A client to the Gene Ontology."""
import logging
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
