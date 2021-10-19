# -*- coding: utf-8 -*-

"""A client to MGI via an ontology client."""

from typing import Optional

from indra.databases.obo_client import PyOboClient

__all__ = [
    "get_name_from_id",
    "get_id_from_name",
]

_client = PyOboClient('mgi')


def get_name_from_id(mgi_id: str) -> Optional[str]:
    """Return the MGI name corresponding to the given MGI.

    Parameters
    ----------
    mgi_id :
        The MGI identifier to be converted. Example: "MGI:97847"

    Returns
    -------
    :
        The MGI name corresponding to the given MGI identifier.

    >>> from indra.databases import mgi_client
    >>> mgi_client.get_name_from_id("MGI:97847")
    'Raf1'
    """
    return _client.get_name_from_id(mgi_id)


def get_id_from_name(mgi_name: str) -> Optional[str]:
    """Return the MGI identifier corresponding to the given MGI name.

    Parameters
    ----------
    mgi_name :
        The MGI name to be converted. Example: "parasite role"

    Returns
    -------
    :
        The MGI identifier corresponding to the given MGI name.

    >>> from indra.databases import mgi_client
    >>> mgi_client.get_name_from_id("Raf1")
    'MGI:97847'
    """
    return _client.get_id_from_name(mgi_name)
