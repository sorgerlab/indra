__all__ = ['standardize_agent_name', 'standardize_db_refs']

import logging
from copy import deepcopy
from indra.ontology.bio import bio_ontology
from indra.statements.agent import default_ns_order

logger = logging.getLogger(__name__)


ns_priorities = {ns: idx for idx, ns in enumerate(default_ns_order)}


def default_prioritize(ns1, ns2):
    ns1p = ns_priorities.get(ns1)
    ns2p = ns_priorities.get(ns2)
    if ns2p is not None and (ns1p is None or ns2p < ns1p):
        return True
    return False


def standardize_db_refs(db_refs, prioritize=default_prioritize):
    """Return a standardized db refs dict for a given db refs dict.

    Parameters
    ----------
    db_refs : dict
        A dict of db refs that may not be standardized, i.e., may be
        missing an available UP ID corresponding to an existing HGNC ID.
    prioritize : function
        A function which takes two str arguments and returns True
        if the second should take priority over the first, and False
        otherwise.

    Returns
    -------
    dict
        The db_refs dict with standardized entries.
    """
    for source_db_ns, source_db_id in deepcopy(db_refs).items():
        mappings = bio_ontology.get_mappings(source_db_ns, source_db_id)
        for mapped_db_ns, mapped_db_id in mappings:
            if mapped_db_ns not in db_refs or \
                    prioritize(mapped_db_ns, source_db_ns):
                db_refs[mapped_db_ns] = mapped_db_id
    return db_refs


def standardize_agent_name(agent, standardize_refs=True):
    """Standardize the name of an Agent based on grounding information.

    If an agent contains a FamPlex grounding, the FamPlex ID is used as a
    name. Otherwise if it contains a Uniprot ID, an attempt is made to find
    the associated HGNC gene name. If one can be found it is used as the
    agent name and the associated HGNC ID is added as an entry to the
    db_refs. Similarly, CHEBI, MESH and GO IDs are used in this order of
    priority to assign a standardized name to the Agent. If no relevant
    IDs are found, the name is not changed.

    Parameters
    ----------
    agent : indra.statements.Agent
        An INDRA Agent whose name attribute should be standardized based
        on grounding information.
    standardize_refs : Optional[bool]
        If True, this function assumes that the Agent's db_refs need to
        be standardized, e.g., HGNC mapped to UP.
        Default: True

    Returns
    -------
    bool
        True if a new name was set, False otherwise.
    """
    # We return immediately for None Agents
    if agent is None:
        return False

    if standardize_refs:
        agent.db_refs = standardize_db_refs(agent.db_refs)

    # We next look for prioritized grounding, if missing, we return
    db_ns, db_id = agent.get_grounding()

    # If there's no grounding then we can't do more to standardize the
    # name and return
    if not db_ns or not db_id:
        return False

    # If there is grounding available, we can try to get the standardized name
    # and in the rare case that we don't get it, we don't set it.
    standard_name = bio_ontology.get_name(db_ns, db_id)
    # Handle special case with UPPRO, if we can't get a feature name
    # we fall back on regular gene/protein naming
    if not standard_name and db_ns == 'UPPRO':
        db_ns, db_id = agent.get_grounding(ns_order=['HGNC', 'UP'])
        if not db_ns or not db_id:
            return False
        standard_name = bio_ontology.get_name(db_ns, db_id)
    if not standard_name:
        return False

    agent.name = standard_name
    return True
