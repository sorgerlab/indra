"""A client to OWL."""

from typing import Optional

from indra.databases.owl_client import OwlClient

_client = OwlClient(prefix='IDO')


def get_ido_name_from_ido_id(ido_id: str) -> Optional[str]:
    """Return the HP name corresponding to the given HP ID.

    Parameters
    ----------
    ido_id : str
        The IDO identifier to be converted. Example: "IDO:0000403"

    Returns
    -------
    :
        The IDO name corresponding to the given IDO identifier.
    """
    return _client.get_name_from_id(ido_id)


def get_ido_id_from_ido_name(ido_name):
    """Return the HP identifier corresponding to the given IDO name.

    Parameters
    ----------
    ido_name : str
        The IDO name to be converted. Example: "parasite role"

    Returns
    -------
    ido_id : str
        The IDO identifier corresponding to the given IDO name.
    """
    return _client.get_id_from_name(ido_name)


if __name__ == '__main__':
    print(*_client.count_xrefs().most_common(), sep='\n')
