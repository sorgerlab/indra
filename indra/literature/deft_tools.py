"""
This file provides several functions helpful for acquiring texts for deft
disambiguation.

It offers the ability to get text content for articles containing a
particular gene. This is useful for aquiring training texts for genes
genes that do not appear in a defining pattern with a problematic shortform.

General XML processing is also provided that allows for extracting
text from a source that may be either of Elsevier XML, NLM XML or raw text.
This is helpful because it avoids having to know in advance the source of
text content from the database.
"""


import time
import logging

from indra.literature import pubmed_client, pmc_client, elsevier_client

logger = logging.getLogger(__name__)

# the elsevier_client will log messages that it is safe to ignore
el = logging.getLogger('indra.literature.elsevier_client')
el.setLevel(logging.WARNING)


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


def get_plaintexts(text_content, contains=None):
    """Returns a corpus of plaintexts given text content from different sources

    Converts xml files into plaintext, leaves abstracts as they are.

    Parameters
    ----------
    sources : list of str
        lists of text content. each item should either be a plaintext, an
        an NLM xml or an Elsevier xml

    Returns
    -------
    plaintexts : list of str
        list of plaintexts for input list of xml strings
    """
    return [universal_extract_text(article, contains)
            for article in text_content]


def universal_extract_text(xml, contains=None):
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
        plaintext = elsevier_client.extract_text(xml, contains)
    except Exception:
        plaintext = None
    if plaintext is None:
        try:
            plaintext = pmc_client.extract_text(xml, contains)
        except Exception:
            plaintext = xml
    return plaintext
