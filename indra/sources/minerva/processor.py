import logging
from indra.ontology.standardize import get_standard_name
from indra.statements import *
from .minerva_client import get_names_to_refs
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
    def __init__(self, sif_strs):
        self.sif_strs = sif_strs
        self.statements = []
        self.names_to_refs = {}

    def extract_statements(self):
        logger.info('Getting names to refs mapping')
        self.names_to_refs = get_names_to_refs()
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
        subj_txt, rel_type, obj_txt = clean_str.split(' ')
        subj = get_agent(subj_txt, self.names_to_refs)
        obj = get_agent(obj_txt, self.names_to_refs)
        if rel_type == 'POSITIVE':
            stmt = Activation(subj, obj)
        elif rel_type == 'NEGATIVE':
            stmt = Inhibition(subj, obj)
        else:
            raise ValueError('Unknown relation type: %s' % rel_type)
        return stmt


def get_agent(txt, names_to_refs):
    txt.replace('_', ' ')
    refs = names_to_refs.get(txt)
    if not refs:
        refs = names_to_refs.get(txt.lower())
    if refs:
        db_refs = indra_db_refs_from_minerva_refs(refs)
        name = get_standard_name(db_refs)
        if not name:
            name = txt
    elif txt.endswith('complex'):
        ag_str = txt[:-8]
        if '/' in ag_str:
            ag_names = ag_str.split('/')
        elif ':' in ag_str:
            ag_names = ag_str.split(':')
        else:
            ag_names = [ag_str]
        agents = [get_agent(ag_name, names_to_refs) for ag_name in ag_names]
        main_agent = agents[0]
        if len(agents) > 1:
            for ag in agents[1:]:
                main_agent.bound_conditions.append(BoundCondition(ag))
        return main_agent
    elif txt.endswith('phenotype'):
        return get_agent(txt[:-10], names_to_refs)
    else:
        name = txt
        db_refs = {'TEXT': txt}
    return Agent(name, db_refs=db_refs)
