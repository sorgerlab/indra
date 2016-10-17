from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
import logging
import itertools
import rdflib.namespace
from copy import deepcopy
from indra.statements import Agent
from indra.preassembler.make_entity_hierarchy import ns_map

logger = logging.getLogger('expand_families')

class Expander(object):
    def __init__(self, hierarchies):
        self.entities = hierarchies['entity']
        # Build reverse lookup dict from the entity hierarchy
        self._children = {}
        logger.info('Generating reverse lookup table for families')
        all_children = set(self.entities.isa_closure.keys()).union(
                       self.entities.partof_closure.keys())
        for child in all_children:
            parents = self._get_parents(child)
            for parent in parents:
                children_list = self._children.get(parent, [])
                children_list.append(child)
                self._children[parent] = children_list

    def _get_parents(self, child):
        immediate_parents = set(self.entities.isa_closure.get(child, [])).union(
                      set(self.entities.partof_closure.get(child, [])))
        all_parents = set([])
        for parent in immediate_parents:
            grandparents = self._get_parents(parent)
            all_parents = all_parents.union(grandparents)
        return all_parents.union(immediate_parents)

    def get_children(self, agent, ns_filter='HGNC'):
        if agent is None:
            return []
        # Get the grounding for the agent
        (ns, id) = agent.get_grounding()
        # Get URI for agent
        ag_uri = self.entities.get_uri(ns, id)
        # Look up the children for this family
        children_uris = self._children.get(ag_uri)
        if children_uris is None:
            return []
        # Parse children URI list into namespaces and ID
        children_parsed = []
        for child_uri in children_uris:
            (child_ns, child_id) = rdflib.namespace.split_uri(child_uri)
            child_ns_name = ns_map.get(child_ns)
            # This shouldn't happen, so generate a logging message if it does
            if child_ns_name is None:
                logger.warning('Unmatched namespace %s in URI %s, skipping' %
                               (child_ns, child_uri))
                continue
            # If ns_filter is None, add in all children
            if ns_filter is None:
                children_parsed.append((child_ns_name, child_id))
            # Otherwise, only add children with a matching namespace
            elif child_ns_name == ns_filter:
                children_parsed.append((child_ns_name, child_id))
        return children_parsed

    def expand_families(self, stmts):
        """Generate statements by expanding members of families and complexes.
        """
        new_stmts = []
        for stmt in stmts:
            # Put together the lists of families, with their members. E.g.,
            # for a statement involving RAF and MEK, should return a list of
            # tuples like [(BRAF, RAF1, ARAF), (MAP2K1, MAP2K2)]
            families_list = []
            for ag in stmt.agent_list():
                ag_children = self.get_children(ag)
                families_list.append(ag_children)
            # Now, put together new statements frmo the cross product of the
            # expanded family members
            for ag_combo in itertools.product(*families_list):
                print(ag_combo)
                # Create agents based on the namespaces/IDs, with appropriate
                # name and db_refs entries
                child_agents = []
                for ag_ns, ag_id in ag_combo:
                    new_agent = Agent(ag_id, db_refs={ag_ns: ag_id})
                    # FIXME FIXME FIXME
                    # This doesn't reproduce agent state from the original
                    # family-level statements!
                    child_agents.append(new_agent)
                # Create a copy of the statement
                new_stmt = deepcopy(stmt)
                # Replace the agents in the statement with the newly-created
                # child agents
                new_stmt.set_agent_list(child_agents)
                print(new_stmt)
                # Add to list
                new_stmts.append(new_stmt)
        return new_stmts
