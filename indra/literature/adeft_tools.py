"""
This file provides several functions helpful for acquiring texts for Adeft
disambiguation.

It offers the ability to get text content for articles containing a
particular gene. This is useful for aquiring training texts for genes
genes that do not appear in a defining pattern with a problematic shortform.

General XML processing is also provided that allows for extracting
text from a source that may be either of Elsevier XML, NLM XML or raw text.
This is helpful because it avoids having to know in advance the source of
text content from the database.
"""

import re
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
            failed.add(pmid)
        time.sleep(0.5)

    remaining_pmids = set(pmids) - pmc_pmids | failed
    abstracts = []
    for pmid in remaining_pmids:
        abstract = pubmed_client.get_abstract(pmid)
        abstracts.append(abstract)
        time.sleep(0.5)

    return [text_content for source in (pmc_xmls, abstracts)
            for text_content in source if text_content is not None]


def universal_extract_paragraphs(xml):
    """Extract paragraphs from xml that could be from  different sources

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
    paragraphs : str
        Extracted plaintext paragraphs from NLM or Elsevier XML
    """
    try:
        paragraphs = elsevier_client.extract_paragraphs(xml)
    except Exception:
        paragraphs = None
    if paragraphs is None:
        try:
            paragraphs = pmc_client.extract_paragraphs(xml)
        except Exception:
            paragraphs = [xml]
    return paragraphs


def filter_paragraphs(paragraphs, contains=None):
    """Filter paragraphs to only those containing one of a list of strings

    Parameters
    ----------
    paragraphs : list of str
        List of plaintext paragraphs from an article

    contains : str or list of str
        Exclude paragraphs not containing this string as a token, or
        at least one of the strings in contains if it is a list

    Returns
    -------
    str
        Plaintext consisting of all input paragraphs containing at least
        one of the supplied tokens.
    """
    if contains is None:
        pattern = ''
    else:
        if isinstance(contains, str):
            contains = [contains]
        pattern = '|'.join(r'[^\w]%s[^\w]' % shortform
                           for shortform in contains)
    paragraphs = [p for p in paragraphs if re.search(pattern, p)]
    return '\n'.join(paragraphs) + '\n'


def universal_extract_text(xml, contains=None):
    """Extract plaintext from xml that could be from different sources

    Parameters
    ----------
    xml : str
        Either an NLM xml, Elsevier xml, or plaintext

    contains : str or list of str
         Exclude paragraphs not containing this string, or at least one
         of the strings in contains if it is a list

    Returns
    -------
    str
        The concatentation of all paragraphs in the input xml, excluding
        paragraphs not containing one of the tokens in the list contains.
        Paragraphs are separated by new lines.
    """
    paragraphs = universal_extract_paragraphs(xml)
    return filter_paragraphs(paragraphs, contains)
