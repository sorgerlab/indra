from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import xml.etree.ElementTree as ET
import os.path
import requests
import logging
from indra.literature import pubmed_client
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('pmc')

pmc_url = 'http://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi'
pmid_convert_url = 'http://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/'

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
    tree = pubmed_client.send_request(pmid_convert_url, data)
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
        logger.warning("Couldn't download PMC%d" % pmc_id)
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
        logger.warning('PMC client returned with error %s: %s' % (err_code, err_text))
        return None
    # If no error, return the XML as a unicode string
    else:
        return xml_bytes.decode('utf-8')


def filter_pmids(pmid_list, source_type):
    """Filter a list of PMIDs for ones with full text from PMC.

    Parameters
    ----------
    pmid_list : list
        List of PMIDs to filter.
    source_type : string
        One of 'fulltext', 'oa_xml', 'oa_txt', or 'auth_xml'.

    Returns
    -------
    list of PMIDs available in the specified source/format type.
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
            fulltext_list = set([line.strip('\n').decode('utf-8')
                                 for line in f.readlines()])
            pmids_fulltext_dict[source_type] = fulltext_list
    return list(set(pmid_list).intersection(
                                pmids_fulltext_dict.get(source_type)))
