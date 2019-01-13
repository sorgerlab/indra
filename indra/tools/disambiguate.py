import time
import logging

from indra.literature.elsevier_client import logger as elsevier_logger
from indra.literature import pubmed_client, pmc_client, elsevier_client

logger = logging.getLogger('disambiguate')

# the elsevier_client will log messages that it is safe to ignore
elsevier_logger.setLevel(logging.WARNING)


def get_text_content_for_gene(hgnc_name):
    """Get articles that have been annotated to contain gene in entrez

    Parameters
    ----------
    hgnc_name : str
       HGNC name for gene

    Returns
    -------
    text_content : list of str
        xmls of fulltext if available otherwise abstracts for all articles
        that haven been annotated in entrez to contain the given gene
    """
    pmids = pubmed_client.get_ids_for_gene(hgnc_name)
    return get_text_content_for_pmids(pmids)


def get_text_content_for_pmids(pmids):
    """Get text content for articles given a list of their pmids

    Parameters
    ----------
    pmids : list of str

    Returns
    -------
    text_content : list of str
    """
    pmc_pmids = set(pmc_client.filter_pmids(pmids, source_type='fulltext'))

    pmc_ids = []
    for pmid in pmc_pmids:
        pmc_id = pmc_client.id_lookup(pmid, idtype='pmid')['pmcid']
        if pmc_id:
            pmc_ids.append(pmc_id)
        else:
            pmc_pmids.discard(pmid)

    pmc_xmls = []
    failed = set()
    for pmc_id in pmc_ids:
        if pmc_id is not None:
            pmc_xmls.append(pmc_client.get_xml(pmc_id))
        else:
            failed.append(pmid)
        time.sleep(0.5)

    remaining_pmids = set(pmids) - pmc_pmids | failed
    abstracts = []
    for pmid in remaining_pmids:
        abstract = pubmed_client.get_abstract(pmid)
        abstracts.append(abstract)
        time.sleep(0.5)

    return [text_content for source in (pmc_xmls, abstracts)
            for text_content in source if text_content is not None]


def build_corpus(*sources):
    """Returns a corpus of plaintexts given text content from different sources

    Converts xml files into plaintext, leaves abstracts as they are.

    Parameters
    ----------
    *sources : list of str
        lists of text content. each item should either be a plaintext, an
        an NLM xml or an Elsevier xml
    """
    return [_universal_extract_text(content) for source in sources
            for content in source]


def _universal_extract_text(xml):
    """Extract plaintext from xml

    First try to parse the xml as if it came from elsevier. if we do not
    have valid elsevier xml this will throw an exception. the text extraction
    function in the pmc client may not throw an exception when parsing elsevier
    xml, silently processing the xml incorrectly

    Parameters
    ----------
    xml : str
       Either an NLM xml, Elsevier xml or plaintext

    Returns
    -------
    plaintext : str
        for NLM or Elsevier xml as input, this is the extracted plaintext
        otherwise the input is returned unchanged
    """
    try:
        plaintext = elsevier_client.extract_text(xml)
    except Exception:
        plaintext = None
    if plaintext is None:
        try:
            plaintext = pmc_client.extract_text(xml)
        except Exception:
            plaintext = xml
    return plaintext
