"""A client to the Disease Ontology."""

from indra.databases.obo_client import OboClient

_client = OboClient(prefix='doid')


def get_doid_name_from_doid_id(doid_id):
    """Return the name corresponding to the given Disease Ontology ID.

    Parameters
    ----------
    doid_id : str
        The Disease Ontology identifier to be converted.
        Example: "DOID:0000017"

    Returns
    -------
    doid_name : str
        The DOID name corresponding to the given DOID identifier.
    """
    return _client.get_name_from_id(doid_id)


def get_doid_id_from_doid_name(doid_name):
    """Return the identifier corresponding to the given Disease Ontology name.

    Parameters
    ----------
    doid_name : str
        The Disease Ontology name to be converted. Example: "Nocturia"

    Returns
    -------
    doid_id : str
        The Disease Ontology identifier corresponding to the given name.
    """
    return _client.get_id_from_name(doid_name)


def get_doid_id_from_doid_alt_id(doid_alt_id):
    """Return the identifier corresponding to the given Disease Ontology alt id.

    Parameters
    ----------
    doid_alt_id : str
        The Disease Ontology alt id to be converted. Example: "DOID:267"

    Returns
    -------
    doid_id : str
        The Disease Ontology identifier corresponding to the given alt id.
    """
    return _client.get_id_from_alt_id(doid_alt_id)
