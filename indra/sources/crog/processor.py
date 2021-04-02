# -*- coding: utf-8 -*-

"""Processor for the `Chemical Roles Graph (CRoG)
<https://github.com/chemical-roles/chemical-roles>`_.
"""

from ..utils import RemoteProcessor

__all__ = [
    'CrogProcessor',
]

CROG_URL = 'https://raw.githubusercontent.com/chemical-roles/' \
           'chemical-roles/master/docs/_data/crog.indra.json'


class CrogProcessor(RemoteProcessor):
    """A processor for INDRA JSON file to be retrieved by URL."""

    def __init__(self):
        super().__init__(url=CROG_URL)
