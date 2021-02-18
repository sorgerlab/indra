import logging
from indra.sources import indra_db_rest
from indra.pipeline import AssemblyPipeline
from indra.sources.hypothesis import upload_statement_annotation

ref_priority = ['TRID', 'PMCID', 'PMID']

logger = logging.getLogger(__name__)


def annotate_paper(text_refs, pipeline=None):
    """Upload INDRA Statements as annotations for a given paper.

    Parameters
    ----------
    text_refs : dict
        A dict of text references, following the same format as
        the INDRA Evidence text_refs attribute.

    pipeline : Optional[json]
        A list of pipeline steps (typically filters) that are applied
        before uploading statements to hypothes.is as annotations.
    """
    for ref_ns in ref_priority:
        ref_id = text_refs.get(ref_ns)
        if ref_id:
            break
    else:
        logger.info('Could not find appropriate text refs')
        return
    ip = indra_db_rest.get_statements_for_paper([(ref_ns.lower(), ref_id)])
    stmts = ip.statements
    # Cut down evidences to ones just from this paper
    for stmt in stmts:
        stmt.evidence = [ev for ev in stmt.evidence if
                         ev.text_refs.get(ref_ns) == ref_id]
    if pipeline:
        ap = AssemblyPipeline(pipeline)
        stmts = ap.run(stmts)

    logger.info('Uploading %d statements to hypothes.is' % len(stmts))
    for stmt in stmts:
        upload_statement_annotation(stmt, annotate_agents=True)