# -*- coding: utf-8 -*-

"""A client to EC-code via an ontology client."""

from typing import List, Optional

from indra.databases.obo_client import PyOboClient

__all__ = [
    "get_name_from_id",
    "get_id_from_name",
]

_client = PyOboClient("ec-code")


def get_name_from_id(ec_code: str) -> Optional[str]:
    """Return the enzyme name corresponding to the given enzyme class code.

    Parameters
    ----------
    ec_code :
        The enzyme class code to be converted. Example: "1.1.1.1"

    Returns
    -------
    :
        The enzyme class name corresponding to the given enzyme class code

    >>> from indra.databases import ec_client
    >>> ec_client.get_name_from_id("1.1.1.1")
    'Alcohol dehydrogenase'
    """
    return _client.get_name_from_id(ec_code)


def get_id_from_name(name: str) -> Optional[str]:
    """Return the enzyme class code corresponding to the given enzyme class name.

    Parameters
    ----------
    name :
        The enzyme name to be converted. Example: "Alcohol dehydrogenase"

    Returns
    -------
    :
        The enzyme class code corresponding to the given enzyme class name.

    >>> from indra.databases import ec_client

    >>> ec_client.get_id_from_name("Alcohol dehydrogenase")
    '1.1.1.1'
    """
    return _client.get_id_from_name(name)


def get_parents(ec_code: str) -> List[str]:
    """Return parents of the given enzyme class code.

    Parameters
    ----------
    ec_code :
        The enzyme class code to looked up. Example: "1.1.1.1"

    Returns
    -------
    :
        The parents of given enzyme class code

    >>> from indra.databases import ec_client
    >>> ec_client.get_parents("1.1.1.1")
    ['1.1.1']
    """
    return _client.get_parents(ec_code)
