import logging
import requests
from indra.sources import indra_db_rest
from indra.literature import pubmed_client
from indra.pipeline import AssemblyPipeline
from indra.statements import stmts_from_json
from indra.sources.hypothesis import upload_statement_annotation

logger = logging.getLogger(__name__)


def annotate_paper_from_db(text_refs, pipeline=None):
    """Upload INDRA Statements as annotations for a given paper based on content
    for that paper in the INDRA DB.

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


def annotate_paper_from_api(text_refs, text_extractor=None, pipeline=None):
    """Read a paper and upload annotations derived from it to hypothes.is.

    Parameters
    ----------
    text_refs : dict
        A dict of text references, following the same format as
        the INDRA Evidence text_refs attribute.
    text_extractor : Optional[function]
        A function which takes the raw content of a website (e.g., HTML)
        and extracts clean text from it to prepare for machine reading.
        This is only used if the text_refs is a URL (e.g., a Wikipedia page),
        it is not used for PMID or PMCID text_refs where content can be
        pre-processed and machine read directly. Default: None
    pipeline : Optional[json]
        A list of pipeline steps (typically filters) that are applied
        before uploading statements to hypothes.is as annotations.
    """
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
        site_content = requests.get(ref_id).text
        if not site_content:
            logger.info('Could not get content from website')
            return
        if text_extractor:
            text = text_extractor(site_content)
        else:
            text = site_content
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
        for ev in stmt.evidence:
            if ref_ns == 'PMID':
                ev.pmid = ref_id
            ev.text_refs[ref_ns] = ref_id
        upload_statement_annotation(stmt, annotate_agents=True)
