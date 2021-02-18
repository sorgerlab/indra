import logging
from indra.sources import indra_db_rest
from indra.pipeline import AssemblyPipeline
from indra.sources.hypothesis import upload_statement_annotation


logger = logging.getLogger(__name__)


def annotate_paper_from_db(text_refs, pipeline=None):
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
    ref_priority = ['TRID', 'PMCID', 'PMID']
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


def annotate_paper_from_api(text_refs, pipeline=None):
    """Read a paper and upload annotations derived from it to hypothes.is.

    Parameters
    ----------
    text_refs : dict
        A dict of text references, following the same format as
        the INDRA Evidence text_refs attribute.

    pipeline : Optional[json]
        A list of pipeline steps (typically filters) that are applied
        before uploading statements to hypothes.is as annotations.
    """
    import requests
    from indra.literature import pubmed_client
    from indra.statements import stmts_from_json
    api_url = 'http://api.indra.bio:8000/reach/'
    ref_priority = ['PMCID', 'PMID', 'URL']
    for ref_ns in ref_priority:
        ref_id = text_refs.get(ref_ns)
        if ref_id:
            break
    else:
        logger.info('Could not find appropriate text refs')
        return
    logger.info('Selected the following paper ID: %s:%s' % (ref_ns, ref_id))
    # Get text content and the read the text
    if ref_ns == 'PMCID':
        res = requests.post(api_url + 'process_pmc', json={'pmc_id': ref_id})
    elif ref_ns == 'PMID':
        abstract = pubmed_client.get_abstract(ref_id)
        if not abstract:
            logger.info('Could not get abstract from PubMed')
            return
        logger.info('Got abstract')
        res = requests.post(api_url + 'process_text', json={'text': abstract})
    elif ref_ns == 'URL':
        text = requests.get(ref_id).text
        if not text:
            logger.info('Could not get text from website')
            return
        res = requests.post(api_url + 'process_text', json={'text': text})
    else:
        return
    jstmts = res.json().get('statements')
    logger.info('Got %d statements from reading' % len(jstmts))
    stmts = stmts_from_json(res.json().get('statements'))

    if pipeline:
        ap = AssemblyPipeline(pipeline)
        stmts = ap.run(stmts)

    logger.info('Uploading %d statements to hypothes.is' % len(stmts))
    for stmt in stmts:
        upload_statement_annotation(stmt, annotate_agents=True)
