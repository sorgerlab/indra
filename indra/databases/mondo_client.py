"""A client to the Monarch Disease Ontology (MONDO)."""

from typing import Optional

from indra.databases.owl_client import OwlClient

# _client = OwlClient(prefix='mondo')


def get_name_from_id(mondo_id: str) -> Optional[str]:
    """Return the name corresponding to the given MONDO ID.

    Parameters
    ----------
    mondo_id : str
        The MONDO identifier to be converted.
        Example: "0002399"

    Returns
    -------
    :
        The MONDO name corresponding to the given MONDO identifier.
    """
    return _client.get_name_from_id(mondo_id)


def get_id_from_name(mondo_name):
    """Return the identifier corresponding to the given MONDO name.

    Parameters
    ----------
    mondo_name : str
        The MONDO name to be converted. Example: "tenosynovial giant cell tumor, localized type"

    Returns
    -------
    :
        The MONDO identifier corresponding to the given name.
    """
    return _client.get_id_from_name(mondo_name)


def get_id_from_alt_id(mondo_alt_id: str) -> Optional[str]:
    """Return the identifier corresponding to the given MONDO alt id.

    Parameters
    ----------
    mondo_alt_id :
        The MONDO alt id to be converted. Example: "0024812"

    Returns
    -------
    :
        The MONDO identifier corresponding to the given alt id.

    >>> from indra.databases import mondo_client
    >>> assert '0024812' == mondo_client.get_id_from_alt_id('0002399')
    """
    return _client.get_id_from_alt_id(mondo_alt_id)


if __name__ == '__main__':
    OwlClient.update_from_obo_library('mondo', remove_prefix=True)
