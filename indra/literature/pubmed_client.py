"""
Search and get metadata for articles in Pubmed.
"""
import csv
import glob
import gzip
import os
import re
import time
import tqdm
import logging
import random
import subprocess
import requests
from time import sleep
from typing import List
from pathlib import Path
from functools import lru_cache
import xml.etree.ElementTree as ET
from indra.resources import RESOURCES_PATH
from indra.util import UnicodeXMLTreeBuilder as UTB
from indra.util import batch_iter, pretty_save_xml


logger = logging.getLogger(__name__)

pubmed_search = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
pubmed_fetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
pubmed_archive = "https://ftp.ncbi.nlm.nih.gov/pubmed"
pubmed_archive_baseline = pubmed_archive + "/baseline/"
pubmed_archive_update = pubmed_archive + "/updatefiles/"
RETRACTIONS_FILE = os.path.join(RESOURCES_PATH, "pubmed_retractions.tsv")


# Send request can't be cached by lru_cache because it takes a dict
# (a mutable/unhashable type) as an argument. We cache the callers instead.
def send_request(url, data, retry_pause=1, max_tries=3):
    try:
        res = requests.get(url, params=data)
    except requests.exceptions.Timeout as e:
        logger.error('PubMed request timed out')
        logger.error('url: %s, data: %s' % (url, data))
        logger.error(e)
        return None
    except requests.exceptions.RequestException as e:
        logger.error('PubMed request exception')
        logger.error('url: %s, data: %s' % (url, data))
        logger.error(e)
        return None
    if res.status_code in {400, 429, 502, 503} and max_tries > 0:
        sleep(retry_pause)
        # Increase the sleep time at random to avoid multiple clients
        # retrying at the same time for e.g. tests
        retry_pause += 0.5 + 1.5 * random.random()
        return send_request(url, data, retry_pause, max_tries - 1)
    if not res.status_code == 200:
        logger.error('Got return code %d from pubmed client.'
                     % res.status_code)
        return None
    tree = ET.XML(res.content, parser=UTB())
    return tree


@lru_cache(maxsize=100)
def get_ids(search_term, **kwargs):
    """Search Pubmed for paper IDs given a search term.

    Search options can be passed as keyword arguments, some of which are
    custom keywords identified by this function, while others are passed on
    as parameters for the request to the PubMed web service
    For details on parameters that can be used in PubMed searches, see
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch Some useful
    parameters to pass are db='pmc' to search PMC instead of pubmed reldate=2
    to search for papers within the last 2 days mindate='2016/03/01',
    maxdate='2016/03/31' to search for papers in March 2016.

    PubMed, by default, limits returned PMIDs to a small number, and this
    number can be controlled by the "retmax" parameter. This function
    uses a retmax value of 10,000 by default (the maximum supported by PubMed)
    that can be changed via the corresponding keyword argument. Note also
    the retstart argument along with retmax to page across batches of IDs.

    PubMed's REST API makes it difficult to retrieve more than 10k
    PMIDs systematically. See the `get_all_ids` function in this module
    for a way to retrieve more than 10k IDs using the PubMed edirect CLI.

    Parameters
    ----------
    search_term : str
        A term for which the PubMed search should be performed.
    use_text_word : Optional[bool]
        If True, the "[tw]" string is appended to the search term to constrain
        the search to "text words", that is words that appear as whole
        in relevant parts of the PubMed entry (excl. for instance the journal
        name or publication date) like the title and abstract. Using this
        option can eliminate spurious search results such as all articles
        published in June for a search for the "JUN" gene, or journal names
        that contain Acad for a search for the "ACAD" gene.
        See also: https://www.nlm.nih.gov/bsd/disted/pubmedtutorial/020_760.html
        Default : True
    kwargs : kwargs
        Additional keyword arguments to pass to the PubMed search as
        parameters.
    """
    use_text_word = kwargs.pop('use_text_word', True)
    if use_text_word:
        search_term += '[tw]'
    params = {'term': search_term,
              'retmax': 10000,
              'retstart': 0,
              'db': 'pubmed',
              'sort': 'pub+date'}
    params.update(kwargs)
    tree = send_request(pubmed_search, params)
    if tree is None:
        return []
    if tree.find('ERROR') is not None:
        logger.error(tree.find('ERROR').text)
        return []
    if tree.find('ErrorList') is not None:
        for err in tree.find('ErrorList'):
            logger.error('Error - %s: %s' % (err.tag, err.text))
        return []
    count = int(tree.find('Count').text)
    id_terms = tree.findall('IdList/Id')
    if id_terms is None:
        return []
    ids = [idt.text for idt in id_terms]
    if count != len(ids):
        logger.warning('Not all ids were retrieved for search %s;\n'
                       'limited at %d.' % (search_term, params['retmax']))
    return ids


def get_id_count(search_term):
    """Get the number of citations in Pubmed for a search query.

    Parameters
    ----------
    search_term : str
        A term for which the PubMed search should be performed.

    Returns
    -------
    int or None
        The number of citations for the query, or None if the query fails.
    """
    params = {'term': search_term,
              'rettype': 'count',
              'db': 'pubmed'}
    tree = send_request(pubmed_search, params)
    if tree is None:
        return None
    else:
        count = list(tree)[0].text
        return int(count)


@lru_cache(maxsize=100)
def get_ids_for_gene(hgnc_name, **kwargs):
    """Get the curated set of articles for a gene in the Entrez database.

    Search parameters for the Gene database query can be passed in as
    keyword arguments.

    Parameters
    ----------
    hgnc_name : str
        The HGNC name of the gene. This is used to obtain the HGNC ID
        (using the hgnc_client module) and in turn used to obtain the Entrez
        ID associated with the gene. Entrez is then queried for that ID.
    """
    from indra.databases import hgnc_client
    # Get the HGNC ID for the HGNC name
    hgnc_id = hgnc_client.get_hgnc_id(hgnc_name)
    if hgnc_id is None:
        raise ValueError('Invalid HGNC name.')
    # Get the Entrez ID
    entrez_id = hgnc_client.get_entrez_id(hgnc_id)
    if entrez_id is None:
        raise ValueError('Entrez ID not found in HGNC table.')
    # Query the Entrez Gene database
    params = {'db': 'gene',
              'retmode': 'xml',
              'id': entrez_id}
    params.update(kwargs)
    tree = send_request(pubmed_fetch, params)
    if tree is None:
        return []
    if tree.find('ERROR') is not None:
        logger.error(tree.find('ERROR').text)
        return []
    # Get all PMIDs from the XML tree
    id_terms = tree.findall('.//PubMedId')
    if id_terms is None:
        return []
    # Use a set to remove duplicate IDs
    ids = list(set([idt.text for idt in id_terms]))
    return ids


def get_ids_for_mesh(mesh_id, major_topic=False, **kwargs):
    """Return PMIDs that are annotated with a given MeSH ID.

    Parameters
    ----------
    mesh_id : str
        The MeSH ID of a term to search for, e.g., D009101.
    major_topic : bool
        If True, only papers for which the given MeSH ID is annotated as
        a major topic are returned. Otherwise all annotations are considered.
        Default: False
    **kwargs
        Any further PudMed search arguments that are passed to
        get_ids.
    """
    from indra.databases import mesh_client
    mesh_name = mesh_client.get_mesh_name(mesh_id)
    if not mesh_name:
        logger.error('Could not get MeSH name for ID %s' % mesh_id)
        return []
    suffix = 'majr' if major_topic else 'mh'
    search_term = '%s [%s]' % (mesh_name, suffix)
    ids = get_ids(search_term, use_text_word=False, **kwargs)
    if mesh_id.startswith('C') and not major_topic:
        # Get pmids for supplementary concepts as well
        search_term = '%s [nm]' % mesh_name
        ids2 = get_ids(search_term, use_text_word=False, **kwargs)
        ids = list(set(ids) | set(ids2))
    return ids


def get_article_xml(pubmed_id):
    """Get the Article subtree a single article from the Pubmed database.

    Parameters
    ----------
    pubmed_id : str
        A PubMed ID.

    Returns
    -------
    xml.etree.ElementTree.Element
        The XML ElementTree Element that represents the Article portion of the
        PubMed entry.
    """
    full_xml_tree = get_full_xml(pubmed_id)
    if full_xml_tree is None:
        return None
    article = full_xml_tree.find('PubmedArticle/MedlineCitation/Article')
    return article  # May be none


@lru_cache(maxsize=100)
def get_full_xml(pubmed_id, fname=None):
    """Get the full XML tree of a single article from the Pubmed database.

    Parameters
    ----------
    pubmed_id : str
        A PubMed ID.
    fname : Optional[str]
        If given, the XML is saved to the given file name.

    Returns
    -------
    xml.etree.ElementTree.Element
        The root element of the XML tree representing the PubMed entry.
        The root is a PubmedArticleSet with a single PubmedArticle element
        that contains the article metadata.
    """
    if pubmed_id.upper().startswith('PMID'):
        pubmed_id = pubmed_id[4:]
    params = {'db': 'pubmed',
              'retmode': 'xml',
              'id': pubmed_id}
    tree = send_request(pubmed_fetch, params)
    if fname:
        pretty_save_xml(tree, fname)
    return tree


def get_title(pubmed_id):
    """Get the title of an article in the Pubmed database."""
    article = get_article_xml(pubmed_id)
    if article is None:
        return None
    return _get_title_from_article_element(article)


def _get_title_from_article_element(article):
    title_tag = article.find('ArticleTitle')
    title = None
    if title_tag is not None:
        title = title_tag.text
        if hasattr(title_tag, 'itertext'):
            title = ''.join(list(title_tag.itertext()))
    return title


def _abstract_from_article_element(article, prepend_title=False):
    abstract = article.findall('Abstract/AbstractText')
    if abstract is None:
        return None
    abstract_text = ' '.join(['' if not hasattr(abst, 'itertext')
                              else ' '.join(list(abst.itertext()))
                              for abst in abstract])
    if prepend_title:
        title = _get_title_from_article_element(article)
        if title is not None:
            if not title.endswith('.'):
                title += '.'
            abstract_text = title + ' ' + abstract_text

    return abstract_text


def get_abstract(pubmed_id, prepend_title=True):
    """Get the abstract of an article in the Pubmed database."""
    article = get_article_xml(pubmed_id)
    if article is None:
        return None
    return _abstract_from_article_element(article, prepend_title)


# A function to get the text for the element, or None if not found
def _find_elem_text(root, xpath_string):
    elem = root.find(xpath_string)
    return None if elem is None else elem.text


def _get_issue_info(journal: ET.Element):
    # Issue info
    issue = journal.find('JournalIssue')
    issue_volume = _find_elem_text(issue, 'Volume')
    issue_issue = _find_elem_text(issue, 'Issue')

    issue_pub_date = issue.find('PubDate')
    if issue_pub_date is not None:
        # Get issue year
        issue_year = _find_elem_text(issue_pub_date, 'Year')
        issue_year = int(issue_year) if issue_year else None

    else:
        issue_year = None

    return {
        "volume": issue_volume,
        "issue": issue_issue,
        "year": issue_year
    }


def get_issn_info(
    medline_citation: ET.Element,
    get_issns_from_nlm: str = "never"
):
    """Given a medline citation, get the issn info from the article

    Parameters
    ----------
    medline_citation : xml.etree.ElementTree.Element
        The MedlineCitation element of the PubMed XML tree.
    get_issns_from_nlm : Literal['never', 'missing', 'always']
        Whether to recover ISSN values from the NLM catalog. Options are
        'never', 'missing', and 'always'. If 'missing', then the ISSN
        values will be recovered from the NLM catalog if they are not found
        in the XML. If 'always', then the ISSN values will be recovered from
        the NLM catalog regardless of whether they are found in the XML.
        Default is 'never' (i.e., never recover from NLM catalog regardless
        of whether they are found in the XML).

    Returns
    -------
    dict
        A dictionary journal, issue, and ISSN info. The structure is as
        follows:
        {
            "journal_title": str,
            "journal_abbrev": str,
            "journal_nlm_id": str,
            "issn_dict": {
                "issn": str,
                "issn_l": str,
                "type": "print"|"electronic"|"other",
            },
            "issue_dict": {
                "volume": str,
                "issue": str,
                "year": int
            }
        }
    """
    if get_issns_from_nlm not in ['never', 'missing', 'always']:
        raise ValueError("get_issns_from_nlm must be one of 'never', "
                         "'missing', or 'always'")
    # Journal info
    journal = medline_citation.find('Article/Journal')
    journal_title = _find_elem_text(journal, 'Title')
    journal_abbrev = _find_elem_text(journal, 'ISOAbbreviation')

    # Issue info
    issue_info = _get_issue_info(journal)

    # Get the ISSN from the article record
    issn_dict = {}
    issn_element = journal.find("ISSN")
    if issn_element is not None:
        issn_type = issn_element.attrib.get("IssnType", "other").lower()
        issn = issn_element.text
        issn_dict["issn"] = issn
        issn_dict["type"] = issn_type

    # Get the linking ISSN from the article record
    issn_linking = _find_elem_text(medline_citation,
                                   "MedlineJournalInfo/ISSNLinking")
    if issn_linking:
        issn_dict["issn_l"] = issn_linking

    nlm_id = _find_elem_text(medline_citation,
                             'MedlineJournalInfo/NlmUniqueID')

    # Get ISSN values from the NLM catalog
    if nlm_id and (
            get_issns_from_nlm == 'always' or
            get_issns_from_nlm == 'missing' and not any(issn_dict.values())
    ):
        nlm_issn_list = get_issns_for_journal(nlm_id)
        if nlm_issn_list:
            issn_dict['alternate_issns'] = nlm_issn_list

    return {
        "journal_title": journal_title,
        "journal_abbrev": journal_abbrev,
        "journal_nlm_id": nlm_id,
        "issn_dict": issn_dict,
        "issue_dict": issue_info,
    }


def _get_journal_info(medline_citation, get_issns_from_nlm: bool):
    # Journal info
    journal = medline_citation.find('Article/Journal')
    journal_title = _find_elem_text(journal, 'Title')
    journal_abbrev = _find_elem_text(journal, 'ISOAbbreviation')

    # Issue info
    issue_info = _get_issue_info(journal)

    # Add the ISSN from the article record
    issn_set = set()
    issn = _find_elem_text(journal, 'ISSN')
    if issn:
        issn_set.add(issn)

    # Add the Linking ISSN from the article record
    issn_linking = _find_elem_text(medline_citation,
                                   'MedlineJournalInfo/ISSNLinking')
    if issn_linking:
        issn_set.add(issn_linking)

    # Now get the list of ISSNs from the NLM Catalog
    nlm_id = _find_elem_text(medline_citation,
                             'MedlineJournalInfo/NlmUniqueID')
    if nlm_id and get_issns_from_nlm:
        nlm_issn_list = get_issns_for_journal(nlm_id)
        if nlm_issn_list:
            issn_set.update(v for _, v in nlm_issn_list)

    # Remove any duplicate issns
    issn_list = list(issn_set)

    return {
        'journal_title': journal_title,
        'journal_abbrev': journal_abbrev,
        'issn_list': issn_list,
        'issn_l': issn_linking,
        'journal_nlm_id': nlm_id,
        'issue': issue_info['issue'],
        'volume': issue_info['volume'],
        'year': issue_info['year'],
    }


def _get_pubmed_publication_date(pubmed_data):
    date_dict = dict.fromkeys(['year', 'month', 'day'])

    # Order potential statuses in order of preferences
    status_list = ['pubmed', 'accepted', 'revised', 'received', 'entrez']

    # Look for various statuses, in order of preference as PubStatus in
    # PubmedPubDate
    for status in status_list:
        pubmed_pub_date = \
                    pubmed_data.find('./History/PubMedPubDate[@PubStatus="%s"]'
                                     % status)
        if pubmed_pub_date is not None:
            break
    else:
        logger.warning("Could not find pub date in: \n%s"
                       % ET.tostring(pubmed_data).decode('utf-8'))
        return date_dict

    def _find_date(element):
        value = _find_elem_text(pubmed_pub_date, element)
        return int(value) if value else None

    # Get date elements from extracted pubmed_pub_date element
    for date_elem in ['Year', 'Month', 'Day']:
        date_dict[date_elem.lower()] = _find_date(date_elem)

    return date_dict


def _parse_author(author_info, include_details=False):
    if not include_details:
        last_name = author_info.find("LastName")
        if last_name is None:
            return None
        return last_name.text

    parsed_info = {
        "last_name": None,
        "first_name": None,
        "initials": None,
        "suffix": None,
        "identifier": None,
        "affiliations": None,
    }
    affiliations = []
    for element in author_info.findall("*"):
        if element.tag == "AffiliationInfo":
            affiliation_name = element.find("Affiliation").text
            identifiers = [e.text for e in element.findall("Identifier")]
            affiliations.append({"name": affiliation_name, "identifiers": identifiers})
        elif element.tag == "LastName":
            parsed_info["last_name"] = element.text
        elif element.tag == "ForeName":
            parsed_info["first_name"] = element.text
        elif element.tag == "Initials":
            parsed_info["initials"] = element.text
        elif element.tag == "Suffix":
            parsed_info["suffix"] = element.text
        elif element.tag == "Identifier":
            parsed_info["identifier"] = element.text
    parsed_info["affiliations"] = affiliations
    return parsed_info


def _get_references(reference_list, only_pmid=True):
    """Return a list of references for an article."""
    if reference_list is None:
        return None

    references = []
    for reference in reference_list.findall('Reference'):
        pmid = _find_elem_text(reference, '*/ArticleId[@IdType="pubmed"]')
        if only_pmid:
            references.append(pmid)
        else:
            ref_dict = {
                'pmid': pmid,
                'doi': _find_elem_text(reference, '*/ArticleId[@IdType="doi"]'),
                'pmcid': _find_elem_text(reference, '*/ArticleId[@IdType="pmcid"]'),
                'citation': _find_elem_text(reference, 'Citation'),
            }
            references.append(ref_dict)
    return references


def _get_article_info(medline_citation, pubmed_data, detailed_authors=False):
    article = medline_citation.find('Article')
    pmid = _find_elem_text(medline_citation, './PMID')
    pii = _find_elem_text(article,
                          './ELocationID[@EIdType="pii"][@ValidYN="Y"]')

    # Look for the DOI in the ELocationID field...
    doi = _find_elem_text(article,
                          './ELocationID[@EIdType="doi"][@ValidYN="Y"]')

    # ...and if that doesn't work, look in the ArticleIdList
    if doi is None:
        doi = _find_elem_text(pubmed_data, './/ArticleId[@IdType="doi"]')

    # Try to get the PMCID
    pmcid = _find_elem_text(pubmed_data, './/ArticleId[@IdType="pmc"]')

    # Title
    title = _get_title_from_article_element(article)

    # Author list
    author_elems = article.findall('AuthorList/Author')
    author_names = None if author_elems is None \
        else [_parse_author(au, detailed_authors) for au in author_elems]

    # Get the page number entry
    page = _find_elem_text(article, 'Pagination/MedlinePgn')

    return {'pmid': pmid, 'pii': pii, 'doi': doi, 'pmcid': pmcid,
            'title': title, 'authors': author_names, 'page': page}


def get_metadata_from_xml_tree(tree, get_issns_from_nlm=False,
                               get_abstracts=False, prepend_title=False,
                               mesh_annotations=True, detailed_authors=False,
                               references_included=None):
    """Get metadata for an XML tree containing PubmedArticle elements.

    Documentation on the XML structure can be found at:
        - https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html
        - https://www.nlm.nih.gov/bsd/licensee/elements_alphabetical.html

    Parameters
    ----------
    tree : xml.etree.ElementTree
        ElementTree containing one or more PubmedArticle elements.
    get_issns_from_nlm : Optional[bool]
        Look up the full list of ISSN number for the journal associated with
        the article, which helps to match articles to CrossRef search results.
        Defaults to False, since it slows down performance.
    get_abstracts : Optional[bool]
        Indicates whether to include the Pubmed abstract in the results.
        Default: False
    prepend_title : Optional[bool]
        If get_abstracts is True, specifies whether the article title should
        be prepended to the abstract text. Default: False
    mesh_annotations : Optional[bool]
        If True, extract mesh annotations from the pubmed entries and include
        in the returned data. If false, don't. Default: True
    detailed_authors : Optional[bool]
        If True, extract as many of the author details as possible, such as
        first name, identifiers, and institutions. If false, only last names
        are returned. Default: False
    references_included : Optional[str]
        If 'detailed', include detailed references in the results. If 'pmid', only include
        the PMID of the reference. If None, don't include references. Default: None

    Returns
    -------
    dict of dicts
        Dictionary indexed by PMID. Each value is a dict containing the
        following fields: 'doi', 'title', 'authors', 'journal_title',
        'journal_abbrev', 'journal_nlm_id', 'issn_list', 'page',
        'volume', 'issue', 'issue_pub_date'.
    """
    # Iterate over the articles and build the results dict
    results = {}
    pm_articles = tree.findall('./PubmedArticle')
    for pm_article in pm_articles:
        medline_citation = pm_article.find('./MedlineCitation')
        pubmed_data = pm_article.find('PubmedData')

        # Build the result
        result = {}
        article_info = _get_article_info(medline_citation, pubmed_data, detailed_authors)
        result.update(article_info)
        journal_info = _get_journal_info(medline_citation, get_issns_from_nlm)
        result.update(journal_info)
        if mesh_annotations:
            context_info = _get_annotations(medline_citation)
            result.update(context_info)
        if references_included:
            references = _get_references(pubmed_data.find('ReferenceList'),
                                         only_pmid=(references_included == 'pmid'))
            result['references'] = references

        publication_date = _get_pubmed_publication_date(pubmed_data)
        result['publication_date'] = publication_date

        # Get the abstracts if requested
        if get_abstracts:
            abstract = _abstract_from_article_element(
                medline_citation.find('Article'),
                prepend_title=prepend_title
                )
            result['abstract'] = abstract

        # Add to dict
        results[article_info['pmid']] = result

    return results


def get_mesh_annotations(pmid):
    """Return a list of MeSH annotations for a given PubMed ID.

    Parameters
    ----------
    pmid : str
        A PubMed ID.

    Returns
    -------
    list of dict
        A list of dicts that represent MeSH annotations with the following keys:
        "mesh" representing the MeSH ID, "text" the standrd name associated with
        the MeSH ID, "major_topic" a boolean flag set depending on whether
        the given MeSH ID is assigned as a major topic to the article, and
        "qualifier" which is a MeSH qualifier ID associated with the annotation,
        if available, otherwise None.
    """
    full_xml_tree = get_full_xml(pmid)
    if not full_xml_tree:
        return None
    medline_citation = full_xml_tree.find('PubmedArticle/MedlineCitation')
    if not medline_citation:
        return None
    annotations = _get_annotations(medline_citation)
    return annotations.get('mesh_annotations')


def _get_annotations(medline_citation):

    def _major_topic(e):
        if e is not None and e.get('MajorTopicYN').upper() == 'Y':
            return True
        return False

    info = []
    for elem in medline_citation.findall('.//MeshHeading'):
        dname = elem.find('DescriptorName')
        qualifier_elems = elem.findall('QualifierName')

        mid = dname.attrib['UI']
        major = _major_topic(dname) or any(_major_topic(qual) for qual
                                           in qualifier_elems)
        qualifiers = [{'text': qual.text, 'mesh': qual.attrib['UI']}
                      for qual in qualifier_elems]
        qual = qualifiers[0] if qualifiers else None

        info.append({'type': 'main', 'mesh': mid, 'text': dname.text,
                     'major_topic': major,
                     # This is only here for backwards compatibility with
                     # INDRA DB which expects a single qualifier or None and
                     # turns the single qualifier into an int internally, so
                     # we can't easily put a joined string of multiple
                     # qualifiers here.
                     'qualifier': qual,
                     # This is the proper full list of qualifiers
                     'qualifiers': qualifiers})
    for elem in medline_citation.findall('.//SupplMeshList/SupplMeshName'):
        info.append({'type': 'supplementary', 'mesh': elem.attrib['UI'], 'text': elem.text,
                     'qualifier': None, 'qualifiers': [],
                     'major_topic': False})
    return {'mesh_annotations': info}


def get_metadata_for_ids(pmid_list, get_issns_from_nlm=False,
                         get_abstracts=False, prepend_title=False,
                         detailed_authors=False, references_included=None):
    """Get article metadata for up to 200 PMIDs from the Pubmed database.

    Parameters
    ----------
    pmid_list : list of str
        Can contain 1-200 PMIDs.
    get_issns_from_nlm : bool
        Look up the full list of ISSN number for the journal associated with
        the article, which helps to match articles to CrossRef search results.
        Defaults to False, since it slows down performance.
    get_abstracts : bool
        Indicates whether to include the Pubmed abstract in the results.
    prepend_title : bool
        If get_abstracts is True, specifies whether the article title should
        be prepended to the abstract text.
    detailed_authors : bool
        If True, extract as many of the author details as possible, such as
        first name, identifiers, and institutions. If false, only last names
        are returned. Default: False
    references_included : Optional[str]
        If 'detailed', include detailed references in the results. If 'pmid', only include
        the PMID of the reference. If None, don't include references. Default: None

    Returns
    -------
    dict of dicts
        Dictionary indexed by PMID. Each value is a dict containing the
        following fields: 'doi', 'title', 'authors', 'journal_title',
        'journal_abbrev', 'journal_nlm_id', 'issn_list', 'page'.
    """
    if len(pmid_list) > 200:
        raise ValueError("Metadata query is limited to 200 PMIDs at a time.")
    params = {'db': 'pubmed',
              'retmode': 'xml',
              'id': pmid_list}
    tree = send_request(pubmed_fetch, params)
    if tree is None:
        return None
    return get_metadata_from_xml_tree(tree, get_issns_from_nlm, get_abstracts,
                                      prepend_title,
                                      detailed_authors=detailed_authors,
                                      references_included=references_included)


def get_metadata_for_all_ids(pmid_list, get_issns_from_nlm=False,
                             get_abstracts=False, prepend_title=False,
                             detailed_authors=False, references_included=None):
    """Get article metadata for up to 200 PMIDs from the Pubmed database.

    Parameters
    ----------
    pmid_list : list of str
        Can contain any number of PMIDs.
    get_issns_from_nlm : bool
        Look up the full list of ISSN number for the journal associated with
        the article, which helps to match articles to CrossRef search results.
        Defaults to False, since it slows down performance.
    get_abstracts : bool
        Indicates whether to include the Pubmed abstract in the results.
    prepend_title : bool
        If get_abstracts is True, specifies whether the article title should
        be prepended to the abstract text.
    detailed_authors : bool
        If True, extract as many of the author details as possible, such as
        first name, identifiers, and institutions. If false, only last names
        are returned. Default: False
    references_included : Optional[str]
        If 'detailed', include detailed references in the results. If 'pmid', only include
        the PMID of the reference. If None, don't include references. Default: None

    Returns
    -------
    dict of dicts
        Dictionary indexed by PMID. Each value is a dict containing the
        following fields: 'doi', 'title', 'authors', 'journal_title',
        'journal_abbrev', 'journal_nlm_id', 'issn_list', 'page'.
    """
    all_metadata = {}
    for ids in tqdm.tqdm(batch_iter(pmid_list, 200), desc='Retrieving metadata'):
        time.sleep(0.1)
        metadata = get_metadata_for_ids(list(ids),
                                        get_issns_from_nlm=get_issns_from_nlm,
                                        get_abstracts=get_abstracts,
                                        prepend_title=prepend_title,
                                        detailed_authors=detailed_authors,
                                        references_included=references_included)
        if metadata is not None:
            all_metadata.update(metadata)
    return all_metadata


@lru_cache(maxsize=1000)
def get_issns_for_journal(nlm_id):
    """Get a dict of the ISSN numbers for a journal given its NLM ID.

    Information on NLM XML DTDs is available at
    https://www.nlm.nih.gov/databases/dtd/
    """
    params = {'db': 'nlmcatalog',
              'retmode': 'xml',
              'id': nlm_id}
    tree = send_request(pubmed_fetch, params)
    if tree is None:
        return None
    issn_list = [(e.attrib.get("IssnType", "other").lower(), e.text)
                 for e in tree.findall('.//ISSN')]
    issn_linking = tree.find('.//ISSNLinking')
    if issn_linking:
        issn_list.append(("linking", issn_linking.text))

    # No ISSNs found!
    if not any(v for k, v in issn_list):
        return None
    return issn_list


def expand_pagination(pages):
    """Convert a page number to long form, e.g., from 456-7 to 456-457."""
    # If there is no hyphen, it's a single page, and we're good to go
    parts = pages.split('-')
    if len(parts) == 1: # No hyphen, so no split
        return pages
    elif len(parts) == 2:
        start = parts[0]
        end = parts[1]
        # If the end is the same number of digits as the start, then we
        # don't change anything!
        if len(start) == len(end):
            return pages
        # Otherwise, replace the last digits of start with the digits of end
        num_end_digits = len(end)
        new_end = start[:-num_end_digits] + end
        return '%s-%s' % (start, new_end)
    else: # More than one hyphen, something weird happened
        logger.warning("Multiple hyphens in page number: %s" % pages)
        return pages


def get_substance_annotations(pubmed_id: str) -> List[str]:
    """Return substance MeSH ID for a given PubMedID.

    Note that substance annotations often overlap with MeSH annotations,
    however, there are cases where a substance annotation is not available
    under MeSH annotations.

    Parameters
    ----------
    pubmed_id :
        PubMedID ID whose substance MeSH ID will be returned.

    Returns
    -------
    :
        Substance MeSH IDs corresponding to the given PubMed paper or
        if None present or a failed query, an empty list will be returned.
    """
    root = get_full_xml(pubmed_id)
    nodes = root.findall('.//MedlineCitation/ChemicalList')
    if len(nodes) == 0:
        logger.error('Could not retrieve substance MeSH IDs for %s' % pubmed_id)
        return []

    uid = [b.attrib.get('UI') for node in nodes
           for c in list(node) for b in c.iter('*')
           if 'UI' in b.attrib]
    return uid


def get_all_ids(search_term):
    """Return all PMIDs for a search term using the edirect CLI.

    This function complements the `get_id` function which uses the PubMed
    REST API but is limited to 10k results and is very difficult to
    generalize to systematically fetch all IDs if there are more than 10k
    results. This function uses the edirect CLI which implements logic
    for paging over results.

    This function only works if edirect is installed and is on your PATH.
    See https://www.ncbi.nlm.nih.gov/books/NBK179288/ for instructions.

    Parameters
    ----------
    search_term : str
        A term for which the PubMed search should be performed.

    Returns
    -------
    list[str]
        A list of PMIDs for the given search term.
    """
    cmd = f'esearch -db pubmed -query "{search_term}" | efetch -format uid'
    res = subprocess.getoutput(cmd)
    # Output is divided by new lines
    elements = res.split('\n')
    # If there are more than 10k IDs, the CLI outputs a . for each
    # iteration, these have to be filtered out
    pmids = [e for e in elements if '.' not in e]
    return pmids


def get_publication_types(article: ET.Element):
    """Return the set of PublicationType for the article

    Parameters
    ----------
    article :
        The XML element for the article. Typically, this is a PubmedArticle
        node.

    Returns
    -------
    : set[str]
        A set of publication type
    """
    return {pt.text for pt in article.find('.//PublicationTypeList')}


def is_retracted(pubmed_id: str) -> bool:
    """Return True if the article with the given PMID has been retracted.

    Parameters
    ----------
    pubmed_id :
        The PMID of the paper to check.

    Returns
    -------
    :
        True if the paper has been retracted, False otherwise.
    """
    return retractions.is_retracted(pubmed_id)


def generate_retractions_file(xml_path: str, download_missing: bool = False):
    """Generate a CSV file of retracted papers from the PubMed XML.

    Parameters
    ----------
    xml_path :
        Path to the directory holding the PubMed XML files. The files will
        be globbed from this directory using the pattern 'pubmed*.xml.gz'.
    download_missing :
        If True, download any missing XML files from the PubMed FTP server.
        Default: False. Note: A full download of the PubMed XML files takes up
        to 5 hours.
    """
    if download_missing:
        ensure_xml_files(xml_path)
    retractions = set()

    files = glob.glob(os.path.join(xml_path, 'pubmed*.xml.gz'))
    if not files:
        raise FileNotFoundError(f"No PubMed XML files found in {xml_path}")

    for xml_file in tqdm.tqdm(files, desc="Processing PubMed XML files"):
        xml_str = gzip.open(xml_file).read()
        tree = ET.XML(xml_str, parser=UTB())
        for article in tree.findall('.//PubmedArticle'):
            pub_types = get_publication_types(article)
            if "Retracted Publication" in pub_types:
                pmid = article.find('.//PMID').text
                retractions.add(pmid)

    if not retractions:
        logger.warning(f"No retractions found from {len(files)} XML files")
        return

    logger.info(f"Writing {len(retractions)} retractions to {RETRACTIONS_FILE}")
    with open(RETRACTIONS_FILE, 'w') as fh:
        fh.write('\n'.join(sorted(retractions)))


def ensure_xml_files(xml_path: str, retries: int = 3):
    """Ensure that the XML files are downloaded and up to date.

    Parameters
    ----------
    xml_path :
        Path to the directory holding the PubMed XML files. The files will
        be globbed from this directory using the pattern 'pubmed*.xml.gz'.
    retries :
        Number of times to retry downloading an individual XML file if there
        is an HTTP error. Default: 3.
    """
    xml_path = Path(xml_path)
    xml_path.mkdir(parents=True, exist_ok=True)

    basefiles = [u for u in _get_urls(pubmed_archive_baseline)]
    updatefiles = [u for u in _get_urls(pubmed_archive_update)]

    # Count successfully downloaded files
    for xml_url in tqdm.tqdm(
            basefiles + updatefiles, desc="Downloading PubMed XML files"
    ):
        xml_file_path = xml_path.joinpath(xml_url.split("/")[-1])
        if not xml_file_path.exists():
            success = _download_xml_gz(xml_url, xml_file_path, retries=retries)
            if not success:
                tqdm.tqdm.write(f"Error downloading {xml_url}, skipping")


def _get_urls(url: str):
    """Get the paths to all XML files on the PubMed FTP server."""
    from bs4 import BeautifulSoup

    logger.info("Getting URL paths from %s" % url)

    # Get page
    response = requests.get(url)
    response.raise_for_status()

    # Make soup
    # Todo: see if it's possible to get the lists of files directly from the
    #  FTP server, rather than scraping the HTML
    soup = BeautifulSoup(response.text, "html.parser")

    # Append trailing slash if not present
    url = url if url.endswith("/") else url + "/"

    # Loop over all links
    for link in soup.find_all("a"):
        href = link.get("href")
        # yield if href matches
        # 'pubmed<2 digit year>n<4 digit file index>.xml.gz'
        # but skip the md5 files
        if href and href.startswith("pubmed") and href.endswith(".xml.gz"):
            yield url + href


def _download_xml_gz(xml_url: str, xml_file: Path, md5_check: bool = True,
                     retries: int = 3) -> bool:
    try:
        resp = requests.get(xml_url)
        resp.raise_for_status()
    except requests.exceptions.RequestException as e:
        if retries > 0:
            tqdm.tqdm.write(f"Error downloading {xml_url}, retrying." + str(e))
            sleep(1)
            return _download_xml_gz(xml_url, xml_file, md5_check, retries - 1)
        else:
            tqdm.tqdm.write(f"Error downloading {xml_url}, skipping")
            return False

    if md5_check:
        from hashlib import md5
        md5_resp = requests.get(xml_url + ".md5")
        checksum = md5(resp.content).hexdigest()
        expected_checksum = re.search(
            r"[0-9a-z]+(?=\n)", md5_resp.content.decode("utf-8")
        ).group()
        if checksum != expected_checksum:
            logger.warning(
                f"Checksum mismatch for {xml_url}, skipping download"
            )
            raise ValueError("Checksum mismatch")

    # Write the file xml.gz file
    with xml_file.open("wb") as fh:
        fh.write(resp.content)

    return True


class Retractions:
    def __init__(self):
        self.retractions = None

    def is_retracted(self, pmid):
        if self.retractions is None:
            with open(RETRACTIONS_FILE, 'r') as fh:
                self.retractions = set(fh.read().splitlines())
        return pmid in self.retractions


retractions = Retractions()
