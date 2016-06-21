import os
import logging
import urllib, urllib2
from functools32 import lru_cache
import xml.etree.ElementTree as ET

logger = logging.getLogger('elsevier')

# THIS FILE IS NOT UNDER VERSION CONTROL
# For more information see http://dev.elsevier.com/
api_key_file = os.path.dirname(os.path.realpath(__file__)) + '/' + \
               'elsevier_api_key'

# Read the API key
try:
    with open(api_key_file, 'rt') as fh:
        api_key = fh.read().strip()
except IOError:
    api_key = None

elsevier_ns = {'dc': 'http://purl.org/dc/elements/1.1/',
               'article': 'http://www.elsevier.com/xml/svapi/article/dtd',
               'ja': 'http://www.elsevier.com/xml/ja/dtd',
               'xocs': 'http://www.elsevier.com/xml/xocs/dtd',
               'common': 'http://www.elsevier.com/xml/common/dtd',
               'atom': 'http://www.w3.org/2005/Atom',
               'prism': 'http://prismstandard.org/namespaces/basic/2.0/'}

@lru_cache(maxsize=100)
def download_article(doi):
    """Download an article in XML format from Elsevier."""
    if doi.lower().startswith('doi:'):
        doi = doi[4:]
    url = 'http://api.elsevier.com/content/article/doi/%s' % doi
    if api_key is None:
        logging.error('Missing API key at %s, could not download article.' %
                      api_key_file)
        return None
    params = {'APIKey': api_key, 'httpAccept': 'text/xml'}
    try:
        res = urllib2.urlopen(url, data=urllib.urlencode(params))
    except urllib2.HTTPError:
        logging.error('Cound not download article %s' % doi)
        return None
    xml = res.read()
    return xml

def get_abstract(doi):
    """Get the abstract of an article from Elsevier."""
    xml = download_article(doi)
    et = ET.fromstring(xml)
    coredata = et.find('article:coredata', elsevier_ns)
    abstract = coredata.find('dc:description', elsevier_ns)
    abs_text = abstract.text
    return abs_text

def get_article(doi, output='txt'):
    """Get the full body of an article from Elsevier. There are two output
    modes: 'txt' strips all xml tags and joins the pieces of text in the main
    text, while 'xml' simply takes the tag containing the body of the article
    and returns it as is . In the latter case, downstream code needs to be
    able to interpret Elsever's XML format. """
    xml = download_article(doi)
    if xml is None:
        return None
    et = ET.fromstring(xml)
    full_text = et.find('article:originalText', elsevier_ns)
    if full_text is None:
        logging.info('Could not find full text for %s.' % doi)
        return None
    main_body = full_text.find('xocs:doc/xocs:serial-item/ja:article/ja:body',
                               elsevier_ns)
    if main_body is None:
        return None
    if output == 'xml':
        return main_body
    elif output == 'txt':
        sections = main_body.findall('common:sections/common:section', elsevier_ns)
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
        logging.error('Unknown output format %s.' % output)
        return None
    return full_txt

@lru_cache(maxsize=100)
def get_dois(query_str, count=100):
    """Search ScienceDirect through the API for articles. See 
    http://api.elsevier.com/content/search/fields/scidir 
    for constructing a query string to pass here.
    Example: 'abstract(BRAF) AND all("colorectal cancer")'
    """
    url = 'http://api.elsevier.com/content/search/scidir'
    if api_key is None:
        logging.error('Missing API key at %s, could not perform search.' %
                      api_key_file)
        return None
    params = {'APIKey': api_key,
              'query': query_str,
              'count': count,
              'httpAccept': 'application/xml',
              'sort': '-coverdate',
              'field': 'doi'}
    res = urllib2.urlopen(url, data=urllib.urlencode(params))
    xml = res.read()
    et = ET.fromstring(xml)
    doi_tags = et.findall('atom:entry/prism:doi', elsevier_ns)
    dois = [dt.text for dt in doi_tags]
    return dois
