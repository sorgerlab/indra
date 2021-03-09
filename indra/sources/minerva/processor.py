import logging
from indra.ontology.standardize import get_standard_name
from indra.ontology.bio import bio_ontology
from indra.statements import *
from .minerva_client import get_ids_to_refs
from .id_mapping import indra_db_refs_from_minerva_refs


logger = logging.getLogger()


class SifProcessor():
    """Processor that extracts INDRA Statements from SIF strings.

    Parameters
    ----------
    model_id_to_sif_strs : dict
        A dictionary mapping a model ID (int) to a list of strings in SIF
        format. Example: {799: ['csa2 POSITIVE sa9', 'csa11 NEGATIVE sa30']}

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements extracted from the SIF strings.
    """
    def __init__(self, model_id_to_sif_strs):
        self.model_id_to_sif_strs = model_id_to_sif_strs
        self.statements = []

    def extract_statements(self):
        for model_id, sif_strs in self.model_id_to_sif_strs.items():
            self.statements += self.process_model(model_id, sif_strs)
        logger.info('Got %d total statements from %d models'
                    % (len(self.statements), len(self.model_id_to_sif_strs)))

    def process_model(self, model_id, sif_strs):
        logger.info('Processing model %d' % model_id)
        ids_to_refs, complex_members = get_ids_to_refs(model_id)
        stmts = []
        for sif_str in sif_strs:
            stmt = self.get_stmt(sif_str, ids_to_refs, complex_members)
            if stmt:
                stmts.append(stmt)
        logger.info('Got %d statements from model %d' % (len(stmts), model_id))
        return stmts

    def get_stmt(self, sif_str, ids_to_refs, complex_members):
        if sif_str.startswith('#') or sif_str == '':
            return
        clean_str = sif_str.strip('\n')
        subj_id, rel_type, obj_id = clean_str.split(' ')
        subj = get_agent(subj_id, ids_to_refs, complex_members)
        obj = get_agent(obj_id, ids_to_refs, complex_members)
        if rel_type == 'POSITIVE':
            stmt = Activation(subj, obj)
        elif rel_type == 'NEGATIVE':
            stmt = Inhibition(subj, obj)
        else:
            raise ValueError('Unknown relation type: %s' % rel_type)
        return stmt


def get_agent(elementId, ids_to_refs, complex_members):
    """Get an agent for a MINERVA element.

    Parameters
    ----------
    elementId : str
        ID of an element used in MINERVA API and raw SIF files.
    ids_to_refs : dict
        A dictionary mapping element IDs to MINERVA provided references. Note
        that this mapping is unique per model (same IDs can be mapped to
        different refs in different models).
    complex_members : dict
        A dictionary mapping element ID of a complex element to element IDs of
        its members.

    Returns
    -------
    agent : indra.statements.agent.Agent
        INDRA agent created from given refs.
    """
    # Get references from MINERVA and filter to accepted namespaces
    accepted_ns = default_ns_order + ['TEXT']
    refs = ids_to_refs.get(elementId)
    filtered_refs = [ref for ref in refs if ref[0] in accepted_ns]
    # If it's a complex and doesn't have complex level grounding
    if elementId in complex_members and len(filtered_refs) == 1:
        # Sort to always have the same main agent
        member_ids = sorted(complex_members[elementId])
        agents = [get_agent(member_id, ids_to_refs, complex_members)
                  for member_id in member_ids]
        # Try to get a FamPlex family
        fam = get_family(agents)
        if fam:
            # Combine TEXT from MINERVA and found FPLX ID
            agent_refs = filtered_refs + [fam]
            return get_agent_from_refs(agent_refs)
        # Otherwise treat a list of agents as an agent with bound conditions
        else:
            main_agent = agents[0]
            if len(agents) > 1:
                for ag in agents[1:]:
                    main_agent.bound_conditions.append(BoundCondition(ag))
            return main_agent
    # Now we have either individual agents or complexes with complex level
    # grounding (e.g. from GO or MESH)
    else:
        return get_agent_from_refs(filtered_refs)


def get_family(agents):
    """Get a FamPlex family if all of its members are given."""
    family_sets = []
    ag_groundings = []
    for ag in agents:
        gr = ag.get_grounding()
        ag_groundings.append(gr)
        parents = bio_ontology.get_parents(*gr)
        families = {p for p in parents if p[0] == 'FPLX'}
        family_sets.append(families)
    common_families = family_sets[0].intersection(*family_sets)
    if not common_families:
        return
    for fam in common_families:
        children = bio_ontology.get_children(*fam)
        # Check if all family members are present
        if set(children) == set(ag_groundings):
            return fam


def get_agent_from_refs(refs):
    """Get an agent given MINERVA refs."""
    db_refs = indra_db_refs_from_minerva_refs(refs)
    name = get_standard_name(db_refs)
    if not name:
        name = db_refs.get('TEXT')
    if name and db_refs:
        return Agent(name, db_refs=db_refs)
