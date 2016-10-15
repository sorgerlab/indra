from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
import logging

logger = logging.getLogger('expand_families')

class Expander(object):
    def __init__(self, hierarchies):
        self.entities = hierarchies['entity']
        # Build reverse lookup dict from the entity hierarchy
        self.families = {}
        self.complexes = {}
        logger.info('Generating reverse lookup table for families')
        for child, parents in self.entities.isa_closure.items():
            for parent in parents:
                children = self.families.get(parent, [])
                children.append(child)
                self.families[parent] = children
        logger.info('Generating reverse lookup table for complexes')
        for subunit, complexes in self.entities.partof_closure.items():
            for complex in complexes:
                subunits = self.complexes.get(complex, [])
                subunits.append(subunit)
                self.complexes[complex] = subunits

    def get_children(self, agent):
        # Get the grounding for the agent
        (ns, id) = agent.get_grounding()
        # Get URI for agent
        ag_uri = self.entities.get_uri(ns, id)
        # Look up the children for this family
        children = self.families.get(ag_uri)
        if children is None:
            return []
        else:
            return children

    def expand_families(self, stmts):
        pass
