from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import logging
import xml.etree.ElementTree as ET
import requests
# Python3
try:
    from functools import lru_cache
# Python2
except ImportError:
    from functools32 import lru_cache
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('elsevier')

# THE API KEY IS NOT UNDER VERSION CONTROL FOR SECURITY
# For more information see http://dev.elsevier.com/
api_key_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'elsevier_api_key')

# Read the API key from the file
try:
    with open(api_key_file, 'rt') as fh:
        api_key = fh.read().strip()
except IOError:
    api_key = None

# THE ELSEVIER API URL: ***MUST BE HTTPS FOR SECURITY***
elsevier_api_url = 'https://api.elsevier.com/content' # <--- HTTPS
elsevier_article_url = '%s/article/doi' % elsevier_api_url
elsevier_search_url = '%s/search/scidir' % elsevier_api_url

# Namespaces for Elsevier XML elements
elsevier_ns = {'dc': 'http://purl.org/dc/elements/1.1/',
               'article': 'http://www.elsevier.com/xml/svapi/article/dtd',
               'ja': 'http://www.elsevier.com/xml/ja/dtd',
               'xocs': 'http://www.elsevier.com/xml/xocs/dtd',
               'common': 'http://www.elsevier.com/xml/common/dtd',
               'atom': 'http://www.w3.org/2005/Atom',
               'prism': 'http://prismstandard.org/namespaces/basic/2.0/'}


def download_article(doi):
    """Download an article in XML format from Elsevier."""
    if doi.lower().startswith('doi:'):
        doi = doi[4:]
    url = '%s/%s' % (elsevier_article_url, doi)
    if api_key is None:
        logger.error('Missing API key at %s, could not download article.' %
                      api_key_file)
        return None
    params = {'APIKey': api_key, 'httpAccept': 'text/xml'}
    res = requests.get(url, params)
    if not res.status_code == 200:
        logger.error('Could not download article %s' % doi)
        return None
    # Parse the content from the stream and then return the tree
    xml_tree = ET.XML(res.content, parser=UTB())
    return xml_tree


def get_abstract(doi):
    """Get the abstract of an article from Elsevier."""
    xml_tree = download_article(doi)
    if xml_tree is None:
        return None
    coredata = xml_tree.find('article:coredata', elsevier_ns)
    abstract = coredata.find('dc:description', elsevier_ns)
    abs_text = abstract.text
    return abs_text


def get_article(doi, output='txt'):
    """Get the full body of an article from Elsevier. There are two output
    modes: 'txt' strips all xml tags and joins the pieces of text in the main
    text, while 'xml' simply takes the tag containing the body of the article
    and returns it as is . In the latter case, downstream code needs to be
    able to interpret Elsever's XML format. """
    xml_tree = download_article(doi)
    if xml_tree is None:
        return None
    full_text = xml_tree.find('article:originalText', elsevier_ns)
    if full_text is None:
        logger.info('Could not find full text for %s.' % doi)
        return None
    main_body = full_text.find('xocs:doc/xocs:serial-item/ja:article/ja:body',
                               elsevier_ns)
    if main_body is None:
        return None
    if output == 'xml':
        return main_body
    elif output == 'txt':
        sections = main_body.findall('common:sections/common:section',
                                     elsevier_ns)
        full_txt = ''
        for s in sections:
            # Paragraphs that are directly under the section
            pars = s.findall('common:para', elsevier_ns)
            # Paragraphs that are under a section within the section
            pars += s.findall('common:section/common:para', elsevier_ns)
            for p in pars:
                # Get the initial string inside the paragraph
                if p.text is not None:
                    full_txt += p.text
                # When there are tags inside the paragraph (for instance
                # references), we need to take those child elements one by one
                # and get the corresponding tail strings and join these. 
                full_txt += ''.join([c.tail if c.tail is not None 
                                     else '' for c in p.getchildren()])
                full_txt += '\n'
    else:
        logger.error('Unknown output format %s.' % output)
        return None
    return full_txt


@lru_cache(maxsize=100)
def get_dois(query_str, count=100):
    """Search ScienceDirect through the API for articles.

    See http://api.elsevier.com/content/search/fields/scidir for constructing a
    query string to pass here.  Example: 'abstract(BRAF) AND all("colorectal
    cancer")'
    """
    url = '%s/%s' % (elsevier_search_url, query_str)
    if api_key is None:
        logger.error('Missing API key at %s, could not perform search.' %
                      api_key_file)
        return None
    params = {'APIKey': api_key,
              'query': query_str,
              'count': count,
              'httpAccept': 'application/xml',
              'sort': '-coverdate',
              'field': 'doi'}
    res = requests.get(url, params)
    if not res.status_code == 200:
        return None
    tree = ET.XML(res.content, parser=UTB())
    doi_tags = tree.findall('atom:entry/prism:doi', elsevier_ns)
    dois = [dt.text for dt in doi_tags]
    return dois
