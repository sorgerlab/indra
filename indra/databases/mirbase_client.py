"""A client to miRBase."""

import os

__all__ = [
    'get_mirbase_id_from_mirbase_name',
    'get_mirbase_name_from_mirbase_id',
    'get_hgnc_id_from_mirbase_id',
    'get_mirbase_id_from_hgnc_id',
    'get_mirbase_id_from_hgnc_symbol',
]

HERE = os.path.dirname(os.path.abspath(__file__))
MIRBASE_FILE = os.path.join(HERE, os.pardir, 'resources', 'mirbase.tsv')


def get_mirbase_name_from_mirbase_id(mirbase_id):
    """Return the miRBase name corresponding to the given miRBase ID.

    Parameters
    ----------
    mirbase_id : str
        The miRBase ID to be converted. Example: "MI0000060"

    Returns
    -------
    mirbase_name : str
        The miRBase name corresponding to the given miRBase ID.
    """
    return _mirbase_id_to_name.get(mirbase_id)


def get_mirbase_id_from_mirbase_name(mirbase_name):
    """Return the miRBase identifier corresponding to the given miRBase name.

    Parameters
    ----------
    mirbase_name : str
        The miRBase ID to be converted. Example: "hsa-mir-19b-2"

    Returns
    -------
    mirbase_id : str
        The miRBase ID corresponding to the given miRBase name.
    """
    return _mirbase_name_to_id.get(mirbase_name)


def get_hgnc_id_from_mirbase_id(mirbase_id):
    """Return the HGNC ID corresponding to the given miRBase ID.

    Parameters
    ----------
    mirbase_id : str
        The miRBase ID to be converted. Example: "MI0000060"

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given miRBase ID.
    """
    return _mirbase_id_to_hgnc_id.get(mirbase_id)


def get_mirbase_id_from_hgnc_id(hgnc_id):
    """Return the HGNC ID corresponding to the given miRBase ID.

    Parameters
    ----------
    hgnc_id : str
        An HGNC identifier to convert to miRBase, if it is indeed
        an miRNA. Example: "31476"

    Returns
    -------
    mirbase_id : str
        The miRBase ID corresponding to the given HGNC ID.
    """
    return _hgnc_id_to_mirbase_id.get(hgnc_id)


def get_mirbase_id_from_hgnc_symbol(hgnc_symbol):
    """Return the HGNC gene symbol corresponding to the given miRBase ID.

    Parameters
    ----------
    hgnc_symbol : str
        An HGNC gene symbol to convert to miRBase, if it is indeed
        an miRNA. Example: "MIR19B2"

    Returns
    -------
    mirbase_id : str
        The miRBase ID corresponding to the given HGNC gene symbol.
    """
    return _hgnc_symbol_to_mirbase_id.get(hgnc_symbol)


def _read():
    """Read the miRBase data into some lookup dictionaries."""
    mirbase_id_to_name = {}
    mirbase_name_to_id = {}
    hgnc_id_to_mirbase_id = {}
    mirbase_id_to_hgnc_id = {}
    hgnc_symbol_to_mirbase_id = {}
    mirbase_id_to_hgnc_symbol = {}

    with open(MIRBASE_FILE) as file:
        next(file)
        for line in file:
            try:
                mirbase_id, mirbase_name, db, identifier, name = \
                                                line.strip().split('\t')
            except ValueError:  # fails on WORMBASE since no names
                continue

            mirbase_id_to_name[mirbase_id] = mirbase_name
            mirbase_name_to_id[mirbase_name] = mirbase_id

            if db == 'HGNC':
                hgnc_id_to_mirbase_id[identifier] = mirbase_id
                mirbase_id_to_hgnc_id[mirbase_id] = identifier
                hgnc_symbol_to_mirbase_id[name] = mirbase_id
                mirbase_id_to_hgnc_symbol[mirbase_id] = name

    return (
        mirbase_id_to_name,
        mirbase_name_to_id,
        hgnc_id_to_mirbase_id,
        mirbase_id_to_hgnc_id,
        hgnc_symbol_to_mirbase_id,
        mirbase_id_to_hgnc_symbol,
    )


(
    _mirbase_id_to_name,
    _mirbase_name_to_id,
    _hgnc_id_to_mirbase_id,
    _mirbase_id_to_hgnc_id,
    _hgnc_symbol_to_mirbase_id,
    _mirbase_id_to_hgnc_symbol,
) = _read()
