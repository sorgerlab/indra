from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import logging
import os.path
import requests
from lxml import etree
from lxml.etree import QName
import xml.etree.ElementTree as ET

from indra.literature import pubmed_client
from indra.util import UnicodeXMLTreeBuilder as UTB

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger(__name__)

pmc_url = 'https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi'
pmid_convert_url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/'

# Paths to resource files
pmids_fulltext_path = os.path.join(os.path.dirname(__file__),
                                   'pmids_fulltext.txt')
pmids_oa_xml_path = os.path.join(os.path.dirname(__file__),
                                 'pmids_oa_xml.txt')
pmids_oa_txt_path = os.path.join(os.path.dirname(__file__),
                                 'pmids_oa_txt.txt')
pmids_auth_xml_path = os.path.join(os.path.dirname(__file__),
                                   'pmids_auth_xml.txt')
# Define global dict containing lists of PMIDs among mineable PMCs
# to be lazily initialized
pmids_fulltext_dict = {}


def id_lookup(paper_id, idtype=None):
    """This function takes a Pubmed ID, Pubmed Central ID, or DOI
    and use the Pubmed ID mapping
    service and looks up all other IDs from one
    of these. The IDs are returned in a dictionary."""
    if idtype is not None and idtype not in ('pmid', 'pmcid', 'doi'):
        raise ValueError("Invalid idtype %s; must be 'pmid', 'pmcid', "
                         "or 'doi'." % idtype)
    if paper_id.upper().startswith('PMC'):
        idtype = 'pmcid'
    # Strip off any prefix
    if paper_id.upper().startswith('PMID'):
        paper_id = paper_id[4:]
    elif paper_id.upper().startswith('DOI'):
        paper_id = paper_id[3:]
    data = {'ids': paper_id}
    if idtype is not None:
        data['idtype'] = idtype
    try:
        tree = pubmed_client.send_request(pmid_convert_url, data)
    except Exception as e:
        logger.error('Error looking up PMID in PMC: %s' % e)
        return {}
    if tree is None:
        return {}
    record = tree.find('record')
    if record is None:
        return {}
    doi = record.attrib.get('doi')
    pmid = record.attrib.get('pmid')
    pmcid = record.attrib.get('pmcid')
    ids = {'doi': doi,
           'pmid': pmid,
           'pmcid': pmcid}
    return ids


def get_ids(search_term, retmax=1000):
    return pubmed_client.get_ids(search_term, retmax=retmax, db='pmc')


def get_xml(pmc_id):
    """Returns XML for the article corresponding to a PMC ID."""
    if pmc_id.upper().startswith('PMC'):
        pmc_id = pmc_id[3:]
    # Request params
    params = {}
    params['verb'] = 'GetRecord'
    params['identifier'] = 'oai:pubmedcentral.nih.gov:%s' % pmc_id
    params['metadataPrefix'] = 'pmc'
    # Submit the request
    res = requests.get(pmc_url, params)
    if not res.status_code == 200:
        logger.warning("Couldn't download %s" % pmc_id)
        return None
    # Read the bytestream
    xml_bytes = res.content
    # Check for any XML errors; xml_str should still be bytes
    tree = ET.XML(xml_bytes, parser=UTB())
    xmlns = "http://www.openarchives.org/OAI/2.0/"
    err_tag = tree.find('{%s}error' % xmlns)
    if err_tag is not None:
        err_code = err_tag.attrib['code']
        err_text = err_tag.text
        logger.warning('PMC client returned with error %s: %s'
                       % (err_code, err_text))
        return None
    # If no error, return the XML as a unicode string
    else:
        return xml_bytes.decode('utf-8')


def extract_text(xml_string):
    """Get plaintext from the body of the given NLM XML string.

    This plaintext consists of all paragraphs returned by
    indra.literature.pmc_client.extract_paragraphs separated
    by newlines and then finally terminated by a newline.
    See the DocString of extract_paragraphs for more information.

    Parameters
    ----------
    xml_string : str
        String containing valid NLM XML.

    Returns
    -------
    str
        Extracted plaintext.
    """
    paragraphs = extract_paragraphs(xml_string)
    if paragraphs:
        return '\n'.join(paragraphs) + '\n'
    else:
        return None


def extract_paragraphs(xml_string):
    """Returns list of paragraphs in an NLM XML.

    This returns a list of the plaintexts for each paragraph and title in
    the input XML, excluding some paragraphs with text that should not
    be relevant to biomedical text processing.

    Relevant text includes titles, abstracts, and the contents of many body
    paragraphs. Within figures, tables, and floating elements, only captions
    are retained (One exception is that all paragraphs within floating
    boxed-text elements are retained. These elements often contain short
    summaries enriched with useful information.) Due to captions, nested
    paragraphs can appear in an NLM XML document. Occasionally there are
    multiple levels of nesting. If nested paragraphs appear in the input
    document their texts are returned in a pre-ordered traversal. The text
    within child paragraphs is not included in the output associated to the
    parent. Each parent appears in the output before its children. All children
    of an element appear before the elements following sibling.

    All tags are removed from each paragraph in the list that is returned.
    LaTeX surrounded by <tex-math> tags is removed entirely.

    Parameters
    ----------
    xml_string : str
        String containing valid NLM XML.

    Returns
    -------
    list of str
        List of extracted paragraphs from the input NLM XML
    """
    output = []
    tree = etree.fromstring(xml_string.encode('utf-8'))
    # Remove namespaces if any exist
    if tree.tag.startswith('{'):
        for element in tree.getiterator():
            element.tag = etree.QName(element).localname
        etree.cleanup_namespaces(tree)
    # Strip out latex
    _remove_elements_by_tag(tree, 'tex-math')
    # Strip out all content in unwanted elements except the captions
    _replace_unwanted_elements_with_their_captions(tree)
    # First process front element. Titles alt-titles and abstracts
    # are pulled from here.
    front_elements = _select_from_top_level(tree, 'front')
    for element in front_elements:
        output.extend(_extract_from_front(element))
    # All paragraphs except those in unwanted elements are extracted
    # from the article body
    body_elements = _select_from_top_level(tree, 'body')
    for element in body_elements:
        output.extend(_extract_from_body(element))
    # Only the body sections of subarticles are processed. All
    # unwanted elements are removed entirely, including captions.
    # Even boxed-text elements are removed.
    subarticles = _select_from_top_level(tree, 'sub-article')
    for element in subarticles:
        output.extend(_extract_from_subarticle(element))
    return output


def filter_pmids(pmid_list, source_type):
    """Filter a list of PMIDs for ones with full text from PMC.

    Parameters
    ----------
    pmid_list : list of str
        List of PMIDs to filter.
    source_type : string
        One of 'fulltext', 'oa_xml', 'oa_txt', or 'auth_xml'.

    Returns
    -------
    list of str
        PMIDs available in the specified source/format type.
    """
    global pmids_fulltext_dict
    # Check args
    if source_type not in ('fulltext', 'oa_xml', 'oa_txt', 'auth_xml'):
        raise ValueError("source_type must be one of: 'fulltext', 'oa_xml', "
                         "'oa_txt', or 'auth_xml'.")
    # Check if we've loaded this type, and lazily initialize
    if pmids_fulltext_dict.get(source_type) is None:
        fulltext_list_path = os.path.join(os.path.dirname(__file__),
                                          'pmids_%s.txt' % source_type)
        with open(fulltext_list_path, 'rb') as f:
            fulltext_list = set([line.strip().decode('utf-8')
                                 for line in f.readlines()])
            pmids_fulltext_dict[source_type] = fulltext_list
    return list(set(pmid_list).intersection(
                                pmids_fulltext_dict.get(source_type)))


def _select_from_top_level(tree, tag):
    """Return top level elements from an NLM XML tree

    Parameters
    ----------
    tree : :py:class:`lxml.etree._Element`
        lxml element for entire tree of a valid NLM XML

    tag : str
        Tag of top level elements to return
    Returns
    -------
    list
        List containing lxml Element objects of selected top level elements
    """
    if tree.tag == 'article':
        article = tree
    else:
        article = tree.xpath('.//article')
        if not len(article):
            raise ValueError('Input XML contains no article element')
        # We make the assumption each NLM XML contains only one article element.
        # If this is not the case, then only the first article will be processed
        article = article[0]
    output = []
    xpath = './%s' % tag
    for element in article.xpath(xpath):
        output.append(element)
    return output


def _extract_from_front(front_element):
    """Return list of titles and paragraphs from front of NLM XML

    Parameters
    ----------
    front_element : :py:class:`lxml.etree._Element`
        etree element for front of a valid NLM XML
    Returns
    -------
    list of str
        List of relevant plain text titles and paragraphs taken from front
        section of NLM XML. These include the article title, alt title,
        and paragraphs within abstracts. Unwanted paragraphs such as
        author statements are excluded.
    """
    output = []
    title_xpath = './article-meta/title-group/article-title'
    alt_title_xpath = './article-meta/title-group/alt-title'
    abstracts_xpath = './article-meta/abstract'
    for element in front_element.xpath(_xpath_union(title_xpath,
                                                    alt_title_xpath,
                                                    abstracts_xpath)):
        if element.tag == 'abstract':
            # Extract paragraphs from abstracts
            output.extend(_extract_paragraphs_from_tree(element))
        else:
            # No paragraphs in titles, Just strip tags
            output.append(' '.join(element.itertext()))
    return output


def _extract_from_body(body_element):
    """Return list of paragraphs from main article body of NLM XML
    """
    return _extract_paragraphs_from_tree(body_element)


def _extract_from_subarticle(subarticle_element):
    """Return list of relevant paragraphs from a subarticle"""
    # Get only body element
    body = subarticle_element.xpath('./body')
    if not body:
        return []
    body = body[0]
    # Remove float elements. From observation these do not appear to
    # contain any meaningful information within sub-articles.
    for element in body.xpath(".//*[@position='float']"):
        element.getparent().remove(element)
    return _extract_paragraphs_from_tree(body)


def _remove_elements_by_tag(tree, *tags):
    """Remove elements with given tags

    Parameters
    ----------
    xml_str : str
        String of valid NLM XML

    Returns
    -------
    str
        Copy of input XML string  with desired elements removed
    """
    bad_xpath = _xpath_union(*['.//%s' % tag for tag in tags])
    for element in tree.xpath(bad_xpath):
        element.getparent().remove(element)


def _replace_unwanted_elements_with_their_captions(tree):
    """Replace an element with its captions"""
    floats_xpath = "//*[@position='float']"
    figs_xpath = './/fig'
    tables_xpath = './/table-wrap'
    unwanted_xpath = _xpath_union(floats_xpath, figs_xpath, tables_xpath)
    unwanted = tree.xpath(unwanted_xpath)
    # Iterating through xpath nodes in reverse leads to processing these
    # nodes from bottom up.
    for element in unwanted[::-1]:
        # Don't remove floats that are boxed-text elements. These often contain
        # useful information
        if element.tag == 'boxed-text':
            continue
        captions = element.xpath('./caption')
        captions_element = etree.Element('captions')
        for caption in captions:
            captions_element.append(caption)
        element.getparent().replace(element, captions_element)


def _retain_only_pars(tree):
    """Strip out all tags except title and p tags"""
    for element in tree.xpath('.//*'):
        if element.tag == 'title':
            element.tag = 'p'
        parent = element.getparent()
        if parent is not None and element.tag != 'p':
            etree.strip_tags(element.getparent(), element.tag)


def _pull_nested_paragraphs_to_top(tree):
    """Flatten nexted paragraphs in pre-ordered traversal"""
    nested_paragraphs = tree.xpath('./p/p')
    while nested_paragraphs:
        last = None
        old_parent = None
        for p in nested_paragraphs:
            parent = p.getparent()
            if parent != old_parent:
                last = parent
            parent.remove(p)
            if p.tail:
                parent.text += ' ' + p.tail
                p.tail = ''
            last.addnext(p)
            last = p
        nested_paragraphs = tree.xpath('./p/p')


def _extract_paragraphs_from_tree(tree):
    # In NLM xml, all plaintext is within <p> tags and <title> tags.
    # There can be formatting tags nested within these tags, but no
    # unwanted elements such as figures and tables appear nested
    # within <p> tags and <title> tags. xpath local-name()= syntax
    # is used to ignore namespaces in the NLM XML. Only elements
    # directly under the input element are included to avoid
    # duplication for nested paragraphs.
    _retain_only_pars(tree)
    _pull_nested_paragraphs_to_top(tree)
    paragraphs = []
    for element in tree.xpath('./p'):
        paragraph = ''.join([x.strip() for x in element.itertext()])
        paragraphs.append(paragraph)
    return paragraphs


def _xpath_union(*xpath_list):
    return ' | '.join(xpath_list)
