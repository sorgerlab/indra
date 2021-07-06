"""A client to EFO."""

from indra.databases.obo_client import OboClient

_client = OboClient(prefix='efo')


def get_efo_name_from_efo_id(efo_id):
    """Return the EFO name corresponding to the given EFO ID.

    Parameters
    ----------
    efo_id : str
        The EFO identifier to be converted. Example: "0005557"

    Returns
    -------
    efo_name : str
        The EFO name corresponding to the given EFO identifier.
    """
    return _client.get_name_from_id(efo_id)


def get_efo_id_from_efo_name(efo_name):
    """Return the EFO identifier corresponding to the given EFO name.

    Parameters
    ----------
    efo_name : str
        The EFO name to be converted. Example: "gum cancer"

    Returns
    -------
    efo_id : str
        The EFO identifier corresponding to the given EFO name.
    """
    return _client.get_id_from_name(efo_name)
