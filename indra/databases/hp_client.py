"""A client to HP."""

from indra.databases.obo_client import OboClient

_client = OboClient(prefix='hp')


def get_hp_name_from_hp_id(hp_id):
    """Return the HP name corresponding to the given HP ID.

    Parameters
    ----------
    hp_id : str
        The HP identifier to be converted. Example: "HP:0000017"

    Returns
    -------
    hp_name : str
        The HP name corresponding to the given HP identifier.
    """
    return _client.get_name_from_id(hp_id)


def get_hp_id_from_hp_name(hp_name):
    """Return the HP identifier corresponding to the given HP name.

    Parameters
    ----------
    hp_name : str
        The HP name to be converted. Example: "Nocturia"

    Returns
    -------
    hp_id : str
        The HP identifier corresponding to the given HP name.
    """
    return _client.get_id_from_name(hp_name)
