__all__ = ['standardize_agent_name', 'standardize_db_refs']

import logging
from copy import deepcopy
from collections import defaultdict
from indra.statements.agent import default_ns_order

logger = logging.getLogger(__name__)


ns_priorities = {ns: idx for idx, ns in enumerate(default_ns_order)}


def default_prioritize(ns1, ns2):
    ns1p = ns_priorities.get(ns1)
    ns2p = ns_priorities.get(ns2)
    if ns2p is not None and (ns1p is None or ns2p < ns1p):
        return True
    return False


def get_mappings_dict(mappings):
    md = defaultdict(list)
    for db_ns, db_id in mappings:
        md[db_ns].append(db_id)
    return md


def standardize_db_refs(db_refs,
                        prioritize=default_prioritize, ontology=None):
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
    ontology : Optional[indra.ontology.IndraOntology]
        An IndraOntology object, if not provided, the default BioOntology
        is used.

    Returns
    -------
    dict
        The db_refs dict with standardized entries.
    """
    if ontology is None:
        from indra.ontology.bio import bio_ontology
        ontology = bio_ontology

    # We iterate over all the db_refs entries that currently exist
    for source_db_ns, source_db_id in deepcopy(db_refs).items():
        # For the entry we get all its xref mappings as a list
        # of tuples and turn it into a dict keyed by namespace
        mappings = get_mappings_dict(
            ontology.get_mappings(source_db_ns, source_db_id))
        # We iterate over these mappings and check if they should
        # be applied
        for mapped_db_ns, mapped_db_ids in mappings.items():
            # If the db_refs doesn't yet contain a mapping for this
            # name space then we always add this mapping. If there
            # is already an entry for this name space then
            # we overwrite it if the source name space is higher
            # priority than the name space being mapped to.
            if mapped_db_ns not in db_refs or \
                    prioritize(mapped_db_ns, source_db_ns):
                db_refs[mapped_db_ns] = sorted(mapped_db_ids)[0]
    return db_refs


def standardize_agent_name(agent, standardize_refs=True, ontology=None):
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
    ontology : Optional[indra.ontology.IndraOntology]
        An IndraOntology object, if not provided, the default BioOntology
        is used.

    Returns
    -------
    bool
        True if a new name was set, False otherwise.
    """
    if ontology is None:
        from indra.ontology.bio import bio_ontology
        ontology = bio_ontology

    # We return immediately for None Agents
    if agent is None:
        return False

    if standardize_refs:
        agent.db_refs = standardize_db_refs(agent.db_refs, ontology=ontology)

    # We next look for prioritized grounding, if missing, we return
    db_ns, db_id = agent.get_grounding()

    # If there's no grounding then we can't do more to standardize the
    # name and return
    if not db_ns or not db_id:
        return False

    # If there is grounding available, we can try to get the standardized name
    # and in the rare case that we don't get it, we don't set it.
    standard_name = ontology.get_name(db_ns, db_id)
    # Handle special case with UPPRO, if we can't get a feature name
    # we fall back on regular gene/protein naming
    if not standard_name and db_ns == 'UPPRO':
        db_ns, db_id = agent.get_grounding(ns_order=['HGNC', 'UP'])
        if not db_ns or not db_id:
            return False
        standard_name = ontology.get_name(db_ns, db_id)
    if not standard_name:
        return False

    agent.name = standard_name
    return True
