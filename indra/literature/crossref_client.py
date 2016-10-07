from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import re
import logging
import requests
from indra.literature import pubmed_client
# Python3
try:
    from functools import lru_cache
# Python2
except ImportError:
    from functools32 import lru_cache


logger = logging.getLogger('crossref')

crossref_url = 'http://api.crossref.org/'
crossref_search_url = 'http://search.crossref.org/dois'

# THIS FILE IS NOT UNDER VERSION CONTROL
# For more information see:
# http://clickthroughsupport.crossref.org/click-through-service-for-researchers/
api_key_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'cr_clickthrough_key')

# Read the API key
api_key = None
try:
    with open(api_key_file, 'rt') as fh:
        api_key = fh.read().strip()
except IOError:
    logger.warning('CrossRef Clickthrough API key could not be found at:')
    logger.warning(api_key_file)
    api_key = None

@lru_cache(maxsize=100)
def get_metadata(doi):
    """Returns the metadata of an article given its DOI from CrossRef
    as a JSON dict"""
    url = crossref_url + 'works/' + doi
    res = requests.get(url)
    if res.status_code != 200:
        logger.warning('Could not get CrossRef metadata, code %d' %
                       res.status_code)
        return None
    raw_message = res.json()
    metadata = raw_message.get('message')
    return metadata

def get_fulltext_links(doi):
    """Return a list of links to the full text of an article given its DOI.
    Each list entry is a dictionary with keys:
    - URL: the URL to the full text
    - content-type: e.g. text/xml or text/plain
    - content-version
    - intended-application: e.g. text-mining
    """
    metadata = get_metadata(doi)
    if metadata is None:
        return None
    links = metadata.get('link')
    return links

def get_publisher(doi):
    metadata = get_metadata(doi)
    if metadata is None:
        return None
    publisher = metadata.get('publisher')
    return publisher

def get_url(doi):
    metadata = get_metadata(doi)
    if metadata is None:
        return None
    url = metadata.get('URL')
    return url

def get_license_links(doi):
    metadata = get_metadata(doi)
    if metadata is None:
        return None
    licenses = metadata.get('license')
    if licenses is None:
        return None
    urls = [l.get('URL') for l in licenses]
    return urls

def doi_query(pmid, search_limit=10):
    """Get the DOI for a PMID by matching CrossRef and Pubmed metadata.

    Searches CrossRef using the article title and then accepts search hits only
    if they have a matching journal ISSN and page number with what is obtained
    from the Pubmed database.
    """
    # Get article metadata from PubMed
    pubmed_meta_dict = pubmed_client.get_metadata_for_ids([pmid],
                                                        get_issns_from_nlm=True)
    if pubmed_meta_dict is None or pubmed_meta_dict.get(pmid) is None:
        logger.warning('No metadata found in Pubmed for PMID%s' % pmid)
        return None
    # The test above ensures we've got this now
    pubmed_meta = pubmed_meta_dict[pmid]
    # Check if we already got a DOI from Pubmed itself!
    if pubmed_meta.get('doi'):
        return pubmed_meta.get('doi')
    # Check for the title, which we'll need for the CrossRef search
    pm_article_title = pubmed_meta.get('title')
    if pm_article_title is None:
        logger.warning('No article title found in Pubmed for PMID%s' % pmid)
        return None
    # Get the ISSN list
    pm_issn_list = pubmed_meta.get('issn_list')
    if not pm_issn_list:
        logger.warning('No ISSNs found in Pubmed for PMID%s' % pmid)
        return None
    # Get the page number
    pm_page = pubmed_meta.get('page')
    if not pm_page:
        logger.warning('No page number found in Pubmed for PMID%s' % pmid)
        return None
    # Now query CrossRef using the title we've got
    url = crossref_search_url
    params = {'q': pm_article_title, 'sort': 'score'}
    res = requests.get(crossref_search_url, params)
    if res.status_code != 200:
        logger.info('Could not get DOI from CrossRef, code %d' % res.status_code)
        return None
    raw_message = res.json()
    mapped_doi = None
    # Iterate over the search results, looking up XREF metadata
    for result_ix, result in enumerate(raw_message):
        if result_ix > search_limit:
            logger.warning('No match found within first %s results, giving up!'
                           % search_limit)
            break
        xref_doi_url = result['doi']
        # Strip the URL prefix off of the DOI
        m = re.match('^http://dx.doi.org/(.*)$', xref_doi_url)
        xref_doi = m.groups()[0]
        # Get the XREF metadata using the DOI
        xref_meta = get_metadata(xref_doi)
        xref_issn_list = xref_meta.get('ISSN')
        xref_page = xref_meta.get('page')
        # If there's no ISSN info for this article, skip to the next result
        if not xref_issn_list:
            logger.warning('No ISSN found for DOI %s, skipping' % xref_doi_url)
            continue
        # If there's no page info for this article, skip to the next result
        if not xref_page:
            logger.warning('No page number found for DOI %s, skipping' %
                          xref_doi_url)
            continue
        # Now check for an ISSN match by looking for the set intersection
        # between the Pubmed ISSN list and the CrossRef ISSN list.
        matching_issns = set(pm_issn_list).intersection(set(xref_issn_list))
        # Before comparing page numbers, regularize the page numbers a bit.
        # Note that we only compare the first page number, since frequently
        # the final page number will simply be missing in one of the data
        # sources. We also canonicalize page numbers of the form '14E' to
        # 'E14' (which is the format used by Pubmed).
        pm_start_page = pm_page.split('-')[0].upper()
        xr_start_page = xref_page.split('-')[0].upper()
        if xr_start_page.endswith('E'):
            xr_start_page = 'E' + xr_start_page[:-1]
        # Now compare the ISSN list and page numbers
        if matching_issns and pm_start_page == xr_start_page:
            # We found a match!
            mapped_doi = xref_doi
            break
        # Otherwise, keep looking through the results...
    # Return a DOI, or None if we didn't find one that met our matching
    # criteria
    return mapped_doi

