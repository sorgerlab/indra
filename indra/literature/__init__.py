from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import requests
import logging
import xml.etree.ElementTree as ET
from indra.literature import pubmed_client
from indra.literature import pmc_client
from indra.literature import crossref_client
from indra.literature import elsevier_client
try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('literature')

def id_lookup(paper_id, idtype):
    """Take an ID of type PMID, PMCID, or DOI and lookup the other IDs.

    If the DOI is not found in Pubmed, try to obtain the DOI by doing a
    reverse-lookup of the DOI in CrossRef using article metadata.

    Parameters
    ----------
    paper_id : string
        ID of the article.
    idtype : 'pmid', 'pmcid', or 'doi
        Type of the ID.
    """
    if idtype not in ('pmid', 'pmcid', 'doi'):
        raise ValueError("Invalid idtype %s; must be 'pmid', 'pmcid', "
                         "or 'doi'." % idtype)

    ids = {'doi': None, 'pmid': None, 'pmcid': None}
    pmc_id_results = pmc_client.id_lookup(paper_id, idtype)
    # Start with the results of the PMC lookup and then override with the
    # provided ID
    ids['pmid'] = pmc_id_results.get('pmid')
    ids['pmcid'] = pmc_id_results.get('pmcid')
    ids['doi'] = pmc_id_results.get('doi')
    ids[idtype] = paper_id
    # If we gave a DOI, then our work is done after looking for PMID and PMCID
    if idtype == 'doi':
        return ids
    # If we gave a PMID or PMCID, we need to check to see if we got a DOI.
    # If we got a DOI back, we're done.
    elif ids.get('doi'):
        return ids
    # If we get here, then we've given PMID or PMCID and don't have a DOI yet.
    # If we gave a PMCID and have neither a PMID nor a DOI, then we'll run
    # into problems later on when we try to the reverse lookup using CrossRef.
    # So we bail here and return what we have (PMCID only) with a warning.
    if ids.get('pmcid') and ids.get('doi') is None and ids.get('pmid') is None:
        logger.warning('%s: PMCID without PMID or DOI' % ids.get('pmcid'))
        return ids
    # To clarify the state of things at this point:
    assert ids.get('pmid') is not None
    assert ids.get('doi') is None
    # As a last result, we try to get the DOI from CrossRef (which internally
    # tries to get the DOI from Pubmed in the process of collecting the
    # necessary metadata for the lookup):
    ids['doi'] = crossref_client.doi_query(ids['pmid'])
    # It may still be None, but at this point there's nothing we can do...
    return ids


def get_full_text(paper_id, idtype, preferred_content_type='text/xml'):
    if preferred_content_type not in \
            ('text/xml', 'text/plain', 'application/pdf'):
        raise ValueError("preferred_content_type must be one of 'text/xml', "
                         "'text/plain', or 'application/pdf'.")
    ids = id_lookup(paper_id, idtype)
    pmcid = ids.get('pmcid')
    pmid = ids.get('pmid')
    doi = ids.get('doi')
    # First try to find paper via PMC
    if pmcid:
        nxml = pmc_client.get_xml(pmcid)
        if nxml:
            return nxml, 'nxml'
    # If we got here, it means we didn't find the full text in PMC, so we'll
    # need either the DOI (for lookup in CrossRef) and/or the PMID (so we
    # can fall back on the abstract. If by some strange turn we have neither,
    # give up now.
    if not doi and not pmid:
        return (None, None)

    # If it does not have PMC NXML then we attempt to obtain the full-text
    # through the CrossRef Click-through API
    if doi:
        # Get publisher
        publisher = crossref_client.get_publisher(doi)

        # First check for whether this is Elsevier--if so, use the Elsevier
        # client directly, because the Clickthrough API key seems unreliable.
        # For now, return as text.
        if publisher == 'Elsevier BV':
            article = elsevier_client.get_article(doi, output='txt')
            if article is not None:
                return (article, 'txt')

        # FIXME FIXME FIXME
        # Because we don't yet have a way to process non-Elsevier content
        # obtained from CrossRef, which includes both XML of unknown format
        # and PDFs, we just comment this section out for now
        # Check if there are any full text links
        links = crossref_client.get_fulltext_links(doi)
        """
        if links:
            headers = {}
            # Set the Cross Ref Clickthrough API key in the header, if we've
            # got one
            if crossref_client.api_key is not None:
                headers['CR-Clickthrough-Client-Token'] = \
                        crossref_client.api_key
            # Utility function to get particular links by content-type
            def lookup_content_type(link_list, content_type):
                content_list = [l.get('URL') for l in link_list
                                if l.get('content-type') == content_type]
                return None if not content_list else content_list[0]
            # First check for what the user asked for
            if lookup_content_type(links, preferred_content_type):
                req = requests.get(lookup_content_type(links,
                                                       preferred_content_type),
                                   headers=headers)
                if req.status_code == 200:
                    req_content_type = req.headers['Content-Type']
                    return req.text, req_content_type
                elif req.status_code == 400:
                    logger.warning('Full text query returned 400 (Bad Request): '
                                  'Perhaps missing CrossRef Clickthrough API '
                                  'key?')
                    return (None, None)
            # Check for XML first
            if lookup_content_type(links, 'text/xml'):
                req = requests.get(lookup_content_type(links, 'text/xml'),
                                   headers=headers)
                if req.status_code == 200:
                    req_content_type = req.headers['Content-Type']
                    return req.text, req_content_type
                elif req.status_code == 400:
                    logger.warning('Full text query returned 400 (Bad Request):'
                                  'Perhaps missing CrossRef Clickthrough API '
                                  'key?')
                    return (None, None)
            # Next, plain text
            elif lookup_content_type(links, 'text/plain'):
                req = requests.get(lookup_content_type(links, 'text/plain'),
                                   headers=headers)
                if req.status_code == 200:
                    req_content_type = req.headers['Content-Type']
                    return req.text, req_content_type
                elif req.status_code == 400:
                    logger.warning('Full text query returned 400 (Bad Request):'
                                  'Perhaps missing CrossRef Clickthrough API '
                                  'key?')
                    return (None, None)
            elif lookup_content_type(links, 'application/pdf'):
                pass
            # Wiley's links are often of content-type 'unspecified'.
            elif lookup_content_type(links, 'unspecified'):
                req = requests.get(lookup_content_type(links, 'unspecified'),
                                   headers=headers)
                if req.status_code == 200:
                    req_content_type = req.headers['Content-Type']
                    return 'foo', req_content_type
                elif req.status_code == 400:
                    logger.warning('Full text query returned 400 (Bad Request):'
                                  'Perhaps missing CrossRef Clickthrough API '
                                  'key?')
                    return (None, None)
                elif req.status_code == 401:
                    logger.warning('Full text query returned 401 (Unauthorized)')
                    return (None, None)
                elif req.status_code == 403:
                    logger.warning('Full text query returned 403 (Forbidden)')
                    return (None, None)
            else:
                raise Exception("Unknown content type(s): %s" % links)
        elif publisher == 'American Society for Biochemistry & Molecular ' \
                          'Biology (ASBMB)':
            url = crossref_client.get_url(doi)
            return get_asbmb_full_text(url)
        """
        # end FIXME FIXME FIXME

        # No full text links and not a publisher we support. We'll have to
        # fall back to the abstract.
        #elif pmid:
        if pmid:
            abstract = pubmed_client.get_abstract(pmid)
            return abstract, 'abstract'

        # We have a useless DOI and no PMID. Give up.
        else:
            return (None, None)
    # We don't have a DOI but we're guaranteed to have a PMID at this point,
    # so we fall back to the abstract:
    else:
        abstract = pubmed_client.get_abstract(pmid)
        return abstract, 'abstract'

    # We'll only get here if we've missed a combination of conditions
    assert False


def get_asbmb_full_text(url):
    # Get the location of the full text PDF from the target URL
    req = requests.get(url)
    if req.status_code != 200:
        logger.warning('ASBMB full text query returned status code %s: URL %s'
                      % (req.status_code, url))
        return (None, None)
    # If we're here that means that we successfully got the paper URL
    xml_str = req.text
    tree = ET.XML(xml_str, parser=UTB())
    fulltext_elem = tree.find('.//{http://www.w3.org/1999/xhtml}meta'
                              '[@name="citation_fulltext_html_url"]')
    # Couldn't find the element containing the full text URL
    if fulltext_elem is None:
        logger.warning("ASBMB full text: couldn't find the full text URL "
                      "element among the meta tags.")
        return (None, None)
    fulltext_url = fulltext_elem.attrib['content']
    # Now, get the full text HTML page
    req2 = requests.get(fulltext_url)
    if req2.status_code != 200:
        logger.warning('ASBMB full text query returned status code %s: URL %s'
                      % (req.status_code, fulltext_url))
        return (None, None)
    # We've got the full text page!
    # Get all the section elements
    xml_str2 = req2.text
    tree2 = ET.XML(xml_str2, parser=UTB())
    return None, None
