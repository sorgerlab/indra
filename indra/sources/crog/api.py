# -*- coding: utf-8 -*-

"""API for the `Chemical Roles Graph (CRoG)
<https://github.com/chemical-roles/chemical-roles>`_.
"""

from typing import Optional

from .processor import CrogProcessor

__all__ = [
    'process_from_web',
]


def process_from_web(url: Optional[str] = None) -> CrogProcessor:
    """Process statements from CRoG over the web.

    Parameters
    ----------
    url :
        An optional URL. If none given, defaults to
        :data:`indra.sources.crog.processor.CROG_URL`.

    Returns
    -------
    processor : CrogProcessor
        A processor with pre-extrated statements.
    """
    processor = CrogProcessor(url=url)
    processor.extract_statements()
    return processor
