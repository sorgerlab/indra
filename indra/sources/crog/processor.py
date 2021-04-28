# -*- coding: utf-8 -*-

"""Processor for the `Chemical Roles Graph (CRoG)
<https://github.com/chemical-roles/chemical-roles>`_.
"""

from typing import Optional

from ..utils import RemoteProcessor

__all__ = [
    'CrogProcessor',
]

CROG_URL = 'https://raw.githubusercontent.com/chemical-roles/' \
           'chemical-roles/master/docs/_data/crog.indra.json'


class CrogProcessor(RemoteProcessor):
    """A processor for the Chemical Roles Graph.

    Parameters
    ----------
    url :
        An optional URL. If none given, defaults to
        :data:`indra.sources.crog.processor.CROG_URL`.
    """

    def __init__(self, url: Optional[str] = None):
        super().__init__(url=url or CROG_URL)

        new_stmts = []
        for stmt in self.statements:
            # We remap the source API to crog to align with the belief model
            for ev in stmt.evidence:
                ev.source_api = 'crog'
            # We also change the name of targets whose names are ECCODEs to
            # have the EC prefix in their name
            for agent in stmt.real_agent_list():
                if agent.name == agent.db_refs.get('ECCODE'):
                    agent.name = 'EC%s' % agent.name
            new_stmts.append(stmt)
        self.statements = new_stmts