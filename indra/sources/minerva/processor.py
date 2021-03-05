import logging
from indra.ontology.standardize import get_standard_name
from indra.statements import *
from .minerva_client import get_ids_to_refs
from .id_mapping import indra_db_refs_from_minerva_refs


logger = logging.getLogger()


class SifProcessor():
    """Processor that extracts INDRA Statements from SIF strings.

    Parameters
    ----------
    sif_strs : list(str)
        A list of strings representing the interactions.
        Example: CASP9 POSITIVE CASP3.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements extracted from the SIF strings.
    names_to_refs : dict
        A dictionary mapping entity names as they are presented in SIF strings
        to their Minerva references.
    """
    def __init__(self, sif_strs, model_id, ids_to_refs=None,
                 complex_members=None):
        self.sif_strs = sif_strs
        self.model_id = model_id
        self.ids_to_refs = ids_to_refs if ids_to_refs else {}
        self.complex_members = complex_members
        self.statements = []

    def extract_statements(self):
        if not self.ids_to_refs or not self.complex_members:
            logger.info('Getting element IDs to refs mapping')
            self.ids_to_refs, self.complex_members = get_ids_to_refs(
                self.model_id)
        logger.info('Getting statements for %d SIF strings'
                    % len(self.sif_strs))
        for sif_str in self.sif_strs:
            stmt = self.get_stmt(sif_str)
            if stmt:
                self.statements.append(stmt)
        logger.info('Got %d statements' % len(self.statements))

    def get_stmt(self, sif_str):
        if sif_str.startswith('#'):
            return
        clean_str = sif_str.strip('\n')
        subj_id, rel_type, obj_id = clean_str.split(' ')
        subj = get_agent(subj_id, self.ids_to_refs, self.complex_members)
        obj = get_agent(obj_id, self.ids_to_refs, self.complex_members)
        if rel_type == 'POSITIVE':
            stmt = Activation(subj, obj)
        elif rel_type == 'NEGATIVE':
            stmt = Inhibition(subj, obj)
        else:
            raise ValueError('Unknown relation type: %s' % rel_type)
        return stmt


def get_agent(elementId, ids_to_refs, complex_members):
    if elementId in ids_to_refs:
        refs = ids_to_refs.get(elementId)
        if refs:
            db_refs = indra_db_refs_from_minerva_refs(refs)
            name = get_standard_name(db_refs)
            if name and db_refs:
                return Agent(name, db_refs=db_refs)
    elif elementId in complex_members:
        member_ids = complex_members[elementId]
        agents = [get_agent(member_id, ids_to_refs, complex_members)
                  for member_id in member_ids]
        main_agent = agents[0]
        if len(agents) > 1:
            for ag in agents[1:]:
                main_agent.bound_conditions.append(BoundCondition(ag))
        return main_agent
