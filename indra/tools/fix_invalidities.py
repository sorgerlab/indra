__all__ = ['fix_invalidities', 'fix_invalidities_db_refs',
           'fix_invalidities_agent', 'fix_invalidities_context',
           'fix_invalidities_stmt', 'fix_invalidities_evidence']

import re
import copy
from typing import List, Mapping
from indra.databases.identifiers import ensure_prefix_if_needed
from indra.statements import Evidence, Statement, Agent, BioContext, \
    Translocation
from indra.pipeline import register_pipeline


@register_pipeline
def fix_invalidities(stmts: List[Statement]) -> List[Statement]:
    """Fix invalidities in a list of Statements.

    Note that in some cases statements can be filtered out if there is a known
    issue to which there is no fix, e.g., a Translocation statements
    missing both location parameters.

    Parameters
    ----------
    stmts :
        A list of INDRA Statements.

    Returns
    -------
    :
        The list of statements with invalidities fixed.
    """
    new_stmts = []
    for stmt in stmts:
        if isinstance(stmt, Translocation) and not stmt.from_location and \
                not stmt.to_location:
            continue
        fix_invalidities_stmt(stmt)
        new_stmts.append(stmt)
    return new_stmts


def fix_invalidities_stmt(stmt: Statement):
    """Fix invalidities of a single INDRA Statement in place."""
    for ev in stmt.evidence:
        fix_invalidities_evidence(ev)
    for agent in stmt.real_agent_list():
        fix_invalidities_agent(agent)


def fix_invalidities_evidence(ev: Evidence):
    """Fix invalidities of a single INDRA Evidence in place."""
    for k, v in copy.deepcopy(ev.text_refs).items():
        if v is None:
            ev.text_refs.pop(k, None)
        elif not k.isupper():
            ev.text_refs.pop(k)
            ev.text_refs[k.upper()] = v

    if ev.pmid and not re.match(r'^\d+$', ev.pmid):
        ev.pmid = None
    if ev.text_refs.get('PMID') and not re.match(r'^\d+$',
                                                 ev.text_refs['PMID']):
        ev.text_refs.pop('PMID', None)

    if ev.pmid is None and ev.text_refs.get('PMID') is not None:
        ev.pmid = ev.text_refs['PMID']
    elif ev.text_refs.get('PMID') is None and ev.pmid is not None:
        ev.text_refs['PMID'] = ev.pmid

    if ev.context is not None:
        fix_invalidities_context(ev.context)


def fix_invalidities_agent(agent: Agent):
    """Fix invalidities of a single INDRA Agent in place."""
    agent.db_refs = fix_invalidities_db_refs(agent.db_refs)


def fix_invalidities_db_refs(db_refs: Mapping[str, str]) -> Mapping[str, str]:
    """Return a fixed version of a db_refs grounding dict."""
    if 'PUBCHEM' in db_refs and \
            db_refs['PUBCHEM'].startswith('CID'):
        db_refs['PUBCHEM'] = \
            db_refs['PUBCHEM'].replace('CID:', '').strip()

    db_refs = {k: v for k, v in db_refs.items()
               if v is not None}

    for k, v in copy.deepcopy(db_refs).items():
        if k == 'CHEMBL' and not v.startswith('CHEMBL'):
            db_refs[k] = 'CHEMBL%s' % v
        elif k == 'ECCODE':
            db_refs['ECCODE'] = db_refs['ECCODE'].replace('.-', '')
        elif k == 'UNIPROT':
            db_refs.pop(k)
            if v.startswith('SL-'):
                db_refs['UPLOC'] = v
            else:
                db_refs['UP'] = v
        elif k == 'UAZ':
            db_refs.pop('UAZ')
            if v.startswith('CVCL'):
                db_refs['CVCL'] = v
        else:
            new_val = ensure_prefix_if_needed(k, v)
            db_refs[k] = new_val
    return db_refs


def fix_invalidities_context(context: BioContext):
    """Fix invalidities of a single INDRA BioContext in place."""
    entries = [context.species, context.cell_line, context.disease,
               context.cell_type, context.organ, context.location]
    for entry in entries:
        if entry is not None:
            entry.db_refs = fix_invalidities_db_refs(entry.db_refs)