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
    """Get text from the body of the given NLM XML string.

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


def extract_titles(xml_string):
    """Get list containing article title and also alt title if one exists

    Parameters
    ----------
    xml_string : str
        String containing valid NLM XML

    Returns
    -------
    list of str
        Singleton list containing article title, or a two element
        list containing the article title and the alt title if
        the latter exists.
    """
    output = []
    tree = etree.fromstring(xml_string.encode('utf-8'))
    title_group_path_elements = ['front', 'article-meta', 'title-group']
    title_group_xpath = '/'
    for tag in title_group_path_elements:
        title_group_xpath += "/*[local-name()='%s']" % tag
    for tag in ['article-title', 'alt-title']:
        title_xpath = title_group_xpath + "/*[local-name()='%s']" % tag
        elements = tree.xpath(title_xpath)
        # Although there should be exaclty one article-title and at most one
        # alt title we handle this as if there could any number of each as
        # a guard against unusual input. The only assumption made is that
        # article-titles appear before alt-titles
        for element in elements:
            output.append(' '.join(element.itertext()))
    return output


def extract_paragraphs(xml_string):
    """Returns list of paragraphs in an NLM XML.

    Parameters
    ----------
    xml_string : str
        String containing valid NLM XML.

    Returns
    -------
    list of str
        List of extracted paragraphs in an NLM XML
    """
    tree = etree.fromstring(xml_string.encode('utf-8'))

    paragraphs = []
    # In NLM xml, all plaintext is within <p> tags and <title> tags.
    # There can be formatting tags nested within these tags, but no
    # unwanted elements such as figures and tables appear nested
    # within <p> tags and <title> tags. xpath local-name()= syntax
    # is used to ignore namespaces in the NLM XML
    for element in tree.xpath("//*[local-name()='p'] | "
                              "//*[local-name()='title']"):
        paragraph = ' '.join(element.itertext())
        paragraphs.append(paragraph)
    return paragraphs


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
