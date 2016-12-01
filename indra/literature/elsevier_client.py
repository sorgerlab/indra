"""
For information on the Elsevier API, see:
  - API Specification: http://dev.elsevier.com/api_docs.html
  - Authentication: https://dev.elsevier.com/tecdoc_api_authentication.html
"""

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
from indra.util import read_unicode_csv
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('elsevier')


# THE ELSEVIER API URL: ***MUST BE HTTPS FOR SECURITY***
elsevier_api_url = 'https://api.elsevier.com/content' # <--- HTTPS
elsevier_article_url = '%s/article/doi' % elsevier_api_url
elsevier_search_url = '%s/search/scidir' % elsevier_api_url
elsevier_entitlement_url = '%s/article/entitlement/doi' % elsevier_api_url

# Namespaces for Elsevier XML elements
elsevier_ns = {'dc': 'http://purl.org/dc/elements/1.1/',
               'article': 'http://www.elsevier.com/xml/svapi/article/dtd',
               'ja': 'http://www.elsevier.com/xml/ja/dtd',
               'xocs': 'http://www.elsevier.com/xml/xocs/dtd',
               'common': 'http://www.elsevier.com/xml/common/dtd',
               'atom': 'http://www.w3.org/2005/Atom',
               'prism': 'http://prismstandard.org/namespaces/basic/2.0/'}

# THE API KEY IS NOT UNDER VERSION CONTROL FOR SECURITY
# For more information see http://dev.elsevier.com/
api_key_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            'elsevier_api_keys')
api_key_env_name = 'ELSEVIER_API_KEY'
inst_key_env_name = 'ELSEVIER_INST_KEY'

# Try to read the API key from a file
try:
    elsevier_keys = dict(read_unicode_csv(api_key_file))
    # Check whether the institution key is present
    if not elsevier_keys.get('X-ELS-Insttoken'):
        logger.info('Optional institution key X-ELS-Insttoken not found in '
                    'elsevier key file.')
    # Check that the API key entry has the right name
    if not elsevier_keys.get('X-ELS-APIKey'):
        logger.error('API key X-ELS-APIKey not found in elsevier key file.')
        elsevier_keys = None
except IOError:
    logger.warning('Elsevier API keys file could not be read, trying '
                   'environment variables $%s and $%s.' %
                   (api_key_env_name, inst_key_env_name))
    logger.debug('Tried key file: %s' % api_key_file)
    # Try the environment variable for the api key. This one is optional,
    # so if it is not found then we just leave it out of the keys dict
    elsevier_keys = {}
    if inst_key_env_name in os.environ:
        elsevier_keys['X-ELS-Insttoken'] = os.environ.get(inst_key_env_name)
    else:
        logger.info('No Elsevier institution key found in environment '
                    'variable %s.' % inst_key_env_name)
    # Try the environment variable for the api key. This one is required, so
    # if it is not found then we set the keys dict to None
    if api_key_env_name in os.environ:
        elsevier_keys['X-ELS-APIKey'] = os.environ.get(api_key_env_name)
    else:
        logger.warning('No Elsevier API key found in environment variable '
                     '%s.' % api_key_env_name)
        elsevier_keys = None


def check_entitlement(doi):
    if elsevier_keys is None:
        logger.error('Missing API key, could not check article entitlement.')
        return False
    if doi.lower().startswith('doi:'):
        doi = doi[4:]
    url = '%s/%s' % (elsevier_entitlement_url, doi)
    params = {'httpAccept': 'text/xml'}
    res = requests.get(url, params, headers=elsevier_keys)
    if not res.status_code == 200:
        logger.error('Could not check entitlements for article %s: '
                     'status code %d' % (doi, res.status_code))
        logger.error('Response content: %s' % res.text)
        return False


def download_article(doi):
    """Download an article in XML format from Elsevier."""
    if elsevier_keys is None:
        logger.error('Missing API key, could not download article.')
        return None
    if doi.lower().startswith('doi:'):
        doi = doi[4:]
    url = '%s/%s' % (elsevier_article_url, doi)
    params = {'httpAccept': 'text/xml'}
    res = requests.get(url, params, headers=elsevier_keys)
    if not res.status_code == 200:
        logger.error('Could not download article %s: status code %d' %
                     (doi, res.status_code))
        logger.error('Elsevier response: %s' % res.text)
        return None
    # Return the XML content as a unicode string, assuming UTF-8 encoding
    return res.content.decode('utf-8')


def get_abstract(doi):
    """Get the abstract of an article from Elsevier."""
    xml_string = download_article(doi)
    if xml_string is None:
        return None
    assert isinstance(xml_string, str)
    # Build XML ElementTree
    xml_tree = ET.XML(xml_string.encode('utf-8'), parser=UTB())
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
    xml_string = download_article(doi)
    text = extract_text(xml_string)
    return text


def extract_text(xml_string):
    if xml_string is None:
        return None
    #with open('/Users/johnbachman/Desktop/elsevier.xml', 'wb') as f:
    #    f.write(xml_string.encode('utf-8'))
    assert isinstance(xml_string, str)
    # Build XML ElementTree
    xml_tree = ET.XML(xml_string.encode('utf-8'), parser=UTB())
    # Look for full text element
    full_text = xml_tree.find('article:originalText', elsevier_ns)
    if full_text is None:
        logger.info('Could not find full text element article:originalText')
        return None
    article_body = _get_article_body(full_text)
    if article_body:
        return article_body
    raw_text = _get_raw_text(full_text)
    if raw_text:
        return raw_text
    #pdf = _get_pdf_attachment(full_text)
    #if pdf:
    #    return pdf
    return None


def _get_pdf_attachment(full_text_elem):
    attachments = full_text_elem.findall('xocs:doc/xocs:meta/'
                                         'xocs:attachment-metadata-doc/'
                                         'xocs:attachments', elsevier_ns)
    for att_elem in attachments:
        web_pdf = att_elem.find('xocs:web-pdf', elsevier_ns)
        if web_pdf is None:
            continue
        # Check for a MAIN pdf
        pdf_purpose = web_pdf.find('xocs:web-pdf-purpose', elsevier_ns)
        if not pdf_purpose.text == 'MAIN':
            continue
        else:
            return 'pdf'
        #locations = web_pdf.findall('xocs:ucs-locator', elsevier_ns)
        #for loc in locations:
        #    logger.info("PDF location: %s" % loc.text)
    #logger.info('Could not find PDF attachment.')
    return None


def _get_article_body(full_text_elem):
    # Look for ja:article
    main_body = full_text_elem.find('xocs:doc/xocs:serial-item/'
                                    'ja:article/ja:body', elsevier_ns)
    if main_body is not None:
        logger.info("Found main body element "
                    "xocs:doc/xocs:serial-item/ja:article/ja:body")
        return _get_sections(main_body)
    logger.info("Could not find main body element "
                "xocs:doc/xocs:serial-item/ja:article/ja:body")
    # If no luck with ja:article, try ja:converted_article
    main_body = full_text_elem.find('xocs:doc/xocs:serial-item/'
                                    'ja:converted-article/ja:body', elsevier_ns)
    if main_body is not None:
        logger.info("Found main body element "
                    "xocs:doc/xocs:serial-item/ja:converted-article/ja:body")
        return _get_sections(main_body)
    logger.info("Could not find main body element "
                "xocs:doc/xocs:serial-item/ja:converted-article/ja:body")
    # If we haven't returned by this point, then return None
    return None


def _get_sections(main_body_elem):
    # Get content sections
    sections = main_body_elem.findall('common:sections/common:section',
                                      elsevier_ns)
    if len(sections) == 0:
        logger.info("Found no sections in main body")
        return None
    # Concatenate the section content
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
    return full_txt


def _get_raw_text(full_text_elem):
    # Look for raw_text
    raw_text = full_text_elem.find('xocs:doc/xocs:rawtext', elsevier_ns)
    if raw_text is None:
        logger.info("Could not find rawtext element xocs:doc/xocs:rawtext")
        return None
    else:
        logger.info("Found rawtext element xocs:doc/xocs:rawtext")
        return raw_text.text


@lru_cache(maxsize=100)
def get_dois(query_str, count=100):
    """Search ScienceDirect through the API for articles.

    See http://api.elsevier.com/content/search/fields/scidir for constructing a
    query string to pass here.  Example: 'abstract(BRAF) AND all("colorectal
    cancer")'
    """
    url = '%s/%s' % (elsevier_search_url, query_str)
    if elsevier_keys is None:
        logger.error('Missing API key at %s, could not perform search.' %
                      api_key_file)
        return None
    params = {'query': query_str,
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
