import logging
import itertools
from copy import deepcopy
from indra.statements import Agent, Complex, Evidence
from indra.ontology.standardize import standardize_agent_name

logger = logging.getLogger(__name__)


class Expander(object):
    def __init__(self, ontology=None):
        if ontology is None:
            from indra.ontology.bio import bio_ontology
            self.ontology = bio_ontology
        else:
            self.ontology = ontology

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
                if ag is None:
                    families_list.append([ag])
                else:
                    db_ns, db_id = ag.get_grounding()
                    if db_ns == 'FPLX':
                        children = self.ontology.get_children(db_ns, db_id)
                        children = [c for c in children if c[0] == 'HGNC']
                        if children:
                            families_list.append(children)
                            continue
                    families_list.append([ag])
            # Now, put together new statements from the cross product of the
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

    def complexes_from_hierarchy(self):
        # Iterate over the partof_closure to determine all of the complexes
        # and all of their members
        fplx_nodes = [n for n in self.ontology.nodes if n[0] == 'FPLX']
        all_complexes = {}
        for fplx_node in fplx_nodes:
            parts = self.ontology.get_ancestors(*fplx_node, {'partof'})
            if not parts:
                continue
            complex_subunits = all_complexes.get(fplx_node, [])
            for part in parts:
                complex_subunits.append(part)
            all_complexes[fplx_node] = complex_subunits

        # Now iterate over all of the complexes and create Complex statements
        complex_stmts = []
        for complex, subunits in all_complexes.items():
            # Create an Evidence object for the statement with the URI of the
            # complex as the source_id
            ev = Evidence(source_api='famplex', source_id=complex)
            subunit_agents = [_agent_from_ns_id(*su) for su in subunits]
            complex_stmt = Complex(subunit_agents, evidence=[ev])
            complex_stmts.append(complex_stmt)
        return complex_stmts

    def expanded_complexes_from_hierarchy(self):
        complex_stmts = self.complexes_from_hierarchy()
        expanded_complexes = self.expand_families(complex_stmts)
        return expanded_complexes


def _agent_from_ns_id(ag_ns, ag_id):
    # Add the ID as a placeholder name
    agent = Agent(ag_id)
    # If we have a proper grounding, add to db_refs
    if ag_id is not None:
        agent.db_refs[ag_ns] = ag_id
    # Now standardize db_refs and set standardized name
    standardize_agent_name(agent, standardize_refs=True)
    agent.db_refs['TEXT'] = agent.name
    return agent
