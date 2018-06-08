from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
import logging
import itertools
import rdflib.namespace
from copy import deepcopy
from indra.preassembler.hierarchy_manager import HierarchyManager, \
    UnknownNamespaceException, hierarchies as default_hierarchies
from indra.databases import hgnc_client
from indra.statements import Agent, Complex, Evidence

logger = logging.getLogger('expand_families')

class Expander(object):
    def __init__(self, hierarchies=None):
        if hierarchies is None:
            self.entities = default_hierarchies['entity']
        else:
            self.entities = hierarchies['entity']

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
                # If the agent has no children, then we use the agent itself
                if len(ag_children) == 0:
                    families_list.append([ag])
                # Otherwise, we add the tuple of namespaces/IDs for the children
                else:
                    families_list.append(ag_children)
            # Now, put together new statements frmo the cross product of the
            # expanded family members
            for ag_combo in itertools.product(*families_list):
                # Create new agents based on the namespaces/IDs, with
                # appropriate name and db_refs entries
                child_agents = []
                for ag_entry in ag_combo:
                    # If we got an agent, or None, that means there were no
                    # children; so we use the original agent rather than
                    # construct a new agent
                    if ag_entry is None or isinstance(ag_entry, Agent):
                        new_agent = ag_entry
                    # Otherwise, create a new agent from the ns/ID
                    elif isinstance(ag_entry, tuple):
                        # FIXME FIXME FIXME
                        # This doesn't reproduce agent state from the original
                        # family-level statements!
                        ag_ns, ag_id = ag_entry
                        new_agent = _agent_from_ns_id(ag_ns, ag_id)
                    else:
                        raise Exception('Unrecognized agent entry type.')
                    # Add agent to our list of child agents
                    child_agents.append(new_agent)
                # Create a copy of the statement
                new_stmt = deepcopy(stmt)
                # Replace the agents in the statement with the newly-created
                # child agents
                new_stmt.set_agent_list(child_agents)
                # Add to list
                new_stmts.append(new_stmt)
        return new_stmts

    def get_children(self, agent, ns_filter='HGNC'):
        if agent is None:
            return []
        # Get the grounding for the agent
        (ns, id) = agent.get_grounding()
        # If there is no grounding for this agent, then return no children
        # (empty list)
        if ns is None or id is None:
            return []
        # Get URI for agent
        ag_uri = self.entities.get_uri(ns, id)
        # Look up the children for this family
        children_uris = self.entities.get_children(ag_uri)
        if not children_uris:
            return []
        # Parse children URI list into namespaces and ID
        children_parsed = []
        for child_uri in children_uris:
            child_ns, child_id = self.entities.ns_id_from_uri(child_uri)
            # If ns_filter is None, add in all children
            if ns_filter is None:
                children_parsed.append((child_ns, child_id))
            # Otherwise, only add children with a matching namespace
            elif child_ns == ns_filter:
                children_parsed.append((child_ns, child_id))
        return children_parsed

    def complexes_from_hierarchy(self):
        # Iterate over the partof_closure to determine all of the complexes
        # and all of their members
        all_complexes = {}
        for subunit, complexes in self.entities.partof_closure.items():
            for complex in complexes:
                complex_subunits = all_complexes.get(complex, [])
                complex_subunits.append(subunit)
                all_complexes[complex] = complex_subunits
        # Now iterate over all of the complexes and create Complex statements
        complex_stmts = []
        for complex, subunits in all_complexes.items():
            # Create an Evidence object for the statement with the URI of the
            # complex as the source_id
            ev = Evidence(source_api='famplex', source_id=complex)
            subunit_agents = [_agent_from_uri(su) for su in subunits]
            complex_stmt = Complex(subunit_agents, evidence=[ev])
            complex_stmts.append(complex_stmt)
        return complex_stmts

    def expanded_complexes_from_hierarchy(self):
        complex_stmts = self.complexes_from_hierarchy()
        expanded_complexes = self.expand_families(complex_stmts)
        return expanded_complexes


def _agent_from_uri(uri):
    ag_ns, ag_id = HierarchyManager.ns_id_from_uri(uri)
    agent = _agent_from_ns_id(ag_ns, ag_id)
    return agent


def _agent_from_ns_id(ag_ns, ag_id):
    ag_name = ag_id
    db_refs = {'TEXT': ag_name}
    if ag_ns == 'HGNC':
        hgnc_id = hgnc_client.get_hgnc_id(ag_id)
        if hgnc_id is not None:
            db_refs['HGNC'] = hgnc_id
            up_id = hgnc_client.get_uniprot_id(hgnc_id)
            if up_id is not None:
                db_refs['UP'] = up_id
    else:
        if ag_id is not None:
            db_refs[ag_ns] = ag_id
    return Agent(ag_name, db_refs=db_refs)
