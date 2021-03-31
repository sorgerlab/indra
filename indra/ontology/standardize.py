__all__ = ['standardize_agent_name', 'standardize_db_refs', 'get_standard_name',
           'standardize_name_db_refs']

import logging
from copy import deepcopy
from collections import defaultdict
from indra.statements.agent import default_ns_order, get_grounding, Agent

logger = logging.getLogger(__name__)


default_ns_priorities = {ns: idx for idx, ns in enumerate(default_ns_order)}


def prioritize(ns1, ns2, ns_order=None):
    ns_priorities = {ns: idx for idx, ns in enumerate(ns_order)} \
        if ns_order is not None else default_ns_priorities
    ns1p = ns_priorities.get(ns1)
    ns2p = ns_priorities.get(ns2)
    if ns2p is not None and (ns1p is None or ns2p < ns1p):
        return True
    return False


def _get_mappings_dict(mappings):
    md = defaultdict(list)
    for db_ns, db_id in mappings:
        md[db_ns].append(db_id)
    return md


def get_standard_agent(name, db_refs, ontology=None, ns_order=None, **kwargs):
    """Get a standard agent based on the name, db_refs, and a any other kwargs.

    name : str
        The name of the agent that may not be standardized.
    db_refs : dict
        A dict of db refs that may not be standardized, i.e., may be
        missing an available UP ID corresponding to an existing HGNC ID.
    ontology : Optional[indra.ontology.IndraOntology]
        An IndraOntology object, if not provided, the default BioOntology
        is used.
    ns_order : Optional[list]
        A list of namespaces which are in order of priority with higher
        priority namespaces appearing earlier in the list.
    kwargs :
        Keyword arguments to pass to :func:`Agent.__init__`.

    Returns
    -------
    Agent
        A standard agent
    """
    ontology = bio_ontology if not ontology else ontology
    standard_name, db_refs = standardize_name_db_refs(db_refs, ontology=ontology, ns_order=ns_order)
    if standard_name:
        name = standard_name
    assert_valid_db_refs(db_refs)
    return Agent(name, db_refs=db_refs, **kwargs)


def standardize_db_refs(db_refs, ontology=None, ns_order=None):
    """Return a standardized db refs dict for a given db refs dict.

    Parameters
    ----------
    db_refs : dict
        A dict of db refs that may not be standardized, i.e., may be
        missing an available UP ID corresponding to an existing HGNC ID.
    ontology : Optional[indra.ontology.IndraOntology]
        An IndraOntology object, if not provided, the default BioOntology
        is used.
    ns_order : Optional[list]
        A list of namespaces which are in order of priority with higher
        priority namespaces appearing earlier in the list.

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
        source_db_id = _preprocess_for_mapping(source_db_ns, source_db_id)
        # For the entry we get all its xref mappings as a list
        # of tuples and turn it into a dict keyed by namespace
        mappings = _get_mappings_dict(
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
                    prioritize(mapped_db_ns, source_db_ns,
                               ns_order=ns_order):
                db_refs[mapped_db_ns] = sorted(mapped_db_ids)[0]
    return db_refs


def _preprocess_for_mapping(db_ns, db_id):
    if db_ns == 'UP' and db_id is not None and '-' in db_id:
        return db_id.split('-')[0]
    return db_id


def standardize_name_db_refs(db_refs, ontology=None, ns_order=None):
    """Return a standardized name and db refs dict for a given db refs dict.

    Parameters
    ----------
    db_refs : dict
        A dict of db refs that may not be standardized, i.e., may be
        missing an available UP ID corresponding to an existing HGNC ID.
    ontology : Optional[indra.ontology.IndraOntology]
        An IndraOntology object, if not provided, the default BioOntology
        is used.
    ns_order : Optional[list]
        A list of namespaces which are in order of priority with higher
        priority namespaces appearing earlier in the list.

    Returns
    -------
    str or None
        The standard name based on the db refs, None if not available.
    dict
        The db_refs dict with standardized entries.
    """
    db_refs = standardize_db_refs(db_refs, ontology=ontology,
                                  ns_order=ns_order)
    name = get_standard_name(db_refs, ontology=ontology, ns_order=ns_order)
    return name, db_refs


def get_standard_name(db_refs, ontology=None, ns_order=None):
    """Return a standardized name for a given db refs dict.

    Parameters
    ----------
    db_refs : dict
        A dict of db refs that may not be standardized, i.e., may be
        missing an available UP ID corresponding to an existing HGNC ID.
    ontology : Optional[indra.ontology.IndraOntology]
        An IndraOntology object, if not provided, the default BioOntology
        is used.
    ns_order : Optional[list]
        A list of namespaces which are in order of priority with higher
        priority namespaces appearing earlier in the list.

    Returns
    -------
    str or None
        The standard name based on the db refs, None if not available.
    """
    if ontology is None:
        from indra.ontology.bio import bio_ontology
        ontology = bio_ontology

    # We next look for prioritized grounding, if missing, we return
    db_ns, db_id = get_grounding(db_refs, ns_order=ns_order)

    # If there's no grounding then we can't do more to standardize the
    # name and return
    if not db_ns or not db_id:
        return None

    # If there is grounding available, we can try to get the standardized name
    # and in the rare case that we don't get it, we don't set it.
    standard_name = ontology.get_name(db_ns, db_id)
    # Handle special case with UPPRO, if we can't get a feature name
    # we fall back on regular gene/protein naming
    if not standard_name and db_ns == 'UPPRO':
        db_ns, db_id = get_grounding(db_refs, ns_order=['HGNC', 'UP'])
        if not db_ns or not db_id:
            return None
        standard_name = ontology.get_name(db_ns, db_id)
    if not standard_name:
        return None

    return standard_name


def standardize_agent_name(agent, standardize_refs=True, ontology=None,
                           ns_order=None):
    """Standardize the name of an Agent based on grounding information.

    The priority of which namespace is used as the bases for the
    standard name depends on

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
    ns_order : Optional[list]
        A list of namespaces which are in order of priority with higher
        priority namespaces appearing earlier in the list.

    Returns
    -------
    bool
        True if a new name was set, False otherwise.
    """
    # If the Agent is None, we return immediately
    if agent is None:
        return False
    # If we want to standardize the Agent's db_refs, we call this now
    if standardize_refs:
        agent.db_refs = standardize_db_refs(agent.db_refs, ontology=ontology)
    # We next try to get a standard name based on the Agent's grounding
    standard_name = get_standard_name(agent.db_refs, ontology=ontology,
                                      ns_order=ns_order)
    # If we got a proper standard name, we apply it
    if standard_name:
        agent.name = standard_name
        return True
    return False
