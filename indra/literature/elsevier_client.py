"""
For information on the Elsevier API, see:
  - API Specification: http://dev.elsevier.com/api_docs.html
  - Authentication: https://dev.elsevier.com/tecdoc_api_authentication.html
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import logging
import textwrap
import datetime
import xml.etree.ElementTree as ET
import requests
from time import sleep
from indra.util import flatten
from indra import has_config, get_config
# Python3
try:
    from functools import lru_cache, wraps
# Python2
except ImportError:
    from functools32 import lru_cache, wraps
from indra.util import read_unicode_csv
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('elsevier')


# THE ELSEVIER API URL: ***MUST BE HTTPS FOR SECURITY***
elsevier_api_url = 'https://api.elsevier.com/content' # <--- HTTPS
elsevier_article_url_fmt = '%s/article/%%s' % elsevier_api_url
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
ELSEVIER_KEYS = None
API_KEY_ENV_NAME = 'ELSEVIER_API_KEY'
INST_KEY_ENV_NAME = 'ELSEVIER_INST_KEY'


def _ensure_api_keys(task_desc, failure_ret=None):
    """Wrap Elsevier methods which directly use the API keys.

    Ensure that the keys are retrieved from the environment or config file when
    first called, and store global scope. Subsequently use globally stashed
    results and check for required ids.
    """
    def check_func_wrapper(func):
        @wraps(func)
        def check_api_keys(*args, **kwargs):
            global ELSEVIER_KEYS
            if ELSEVIER_KEYS is None:
                ELSEVIER_KEYS = {}
                # Try to read in Elsevier API keys. For each key, first check
                # the environment variables, then check the INDRA config file.
                if not has_config(INST_KEY_ENV_NAME):
                    logger.warning('Institution API key %s not found in config '
                                   'file or environment variable: this will '
                                   'limit access for %s'
                                   % (INST_KEY_ENV_NAME, task_desc))
                ELSEVIER_KEYS['X-ELS-Insttoken'] = get_config(INST_KEY_ENV_NAME)

                if not has_config(API_KEY_ENV_NAME):
                    logger.error('API key %s not found in configuration file '
                                 'or environment variable: cannot %s'
                                 % (API_KEY_ENV_NAME, task_desc))
                    return failure_ret
                ELSEVIER_KEYS['X-ELS-APIKey'] = get_config(API_KEY_ENV_NAME)
            elif 'X-ELS-APIKey' not in ELSEVIER_KEYS.keys():
                logger.error('No Elsevier API key %s found: cannot %s'
                             % (API_KEY_ENV_NAME, task_desc))
                return failure_ret
            return func(*args, **kwargs)
        return check_api_keys
    return check_func_wrapper


@_ensure_api_keys('check article entitlement', False)
def check_entitlement(doi):
    """Check whether IP and credentials enable access to content for a doi.

    This function uses the entitlement endpoint of the Elsevier API to check
    whether an article is available to a given institution. Note that this
    feature of the API is itself not available for all institution keys.
    """
    if doi.lower().startswith('doi:'):
        doi = doi[4:]
    url = '%s/%s' % (elsevier_entitlement_url, doi)
    params = {'httpAccept': 'text/xml'}
    res = requests.get(url, params, headers=ELSEVIER_KEYS)
    if not res.status_code == 200:
        logger.error('Could not check entitlements for article %s: '
                     'status code %d' % (doi, res.status_code))
        logger.error('Response content: %s' % res.text)
        return False
    return True


@_ensure_api_keys('download article')
def download_article(id_val, id_type='doi', on_retry=False):
    """Low level function to get an XML article for a particular id.

    Parameters
    ----------
    id_val : str
        The value of the id.
    id_type : str
        The type of id, such as pmid (a.k.a. pubmed_id), doi, or eid.
    on_retry : bool
        This function has a recursive retry feature, and this is the only time
        this parameter should be used.

    Returns
    -------
    content : str or None
        If found, the content string is returned, otherwise, None is returned.
    """
    if id_type == 'pmid':
        id_type = 'pubmed_id'
    url = '%s/%s' % (elsevier_article_url_fmt % id_type, id_val)
    params = {'httpAccept': 'text/xml'}
    res = requests.get(url, params, headers=ELSEVIER_KEYS)
    if res.status_code == 404:
        logger.info("Resource for %s not available on elsevier." % url)
        return None
    elif res.status_code == 429:
        if not on_retry:
            logger.warning("Broke the speed limit. Waiting half a second then "
                           "trying again...")
            sleep(0.5)
            return download_article(id_val, id_type, True)
        else:
            logger.error("Still breaking speed limit after waiting.")
            logger.error("Elsevier response: %s" % res.text)
            return None
    elif res.status_code != 200:
        logger.error('Could not download article %s: status code %d' %
                     (url, res.status_code))
        logger.error('Elsevier response: %s' % res.text)
        return None
    else:
        content_str = res.content.decode('utf-8')
        if content_str.startswith('<service-error>'):
            logger.error('Got a service error with 200 status: %s'
                         % content_str)
            return None
    # Return the XML content as a unicode string, assuming UTF-8 encoding
    return content_str


def download_article_from_ids(**id_dict):
    """Download an article in XML format from Elsevier matching the set of ids.

    Parameters
    ----------
    <id_type> : str
        You can enter any combination of eid, doi, pmid, and/or pii. Ids will be
        checked in that order, until either content has been found or all ids
        have been checked.

    Returns
    -------
    content : str or None
        If found, the content is returned as a string, otherwise None is
        returned.
    """
    valid_id_types = ['eid', 'doi', 'pmid', 'pii']
    assert all([k in valid_id_types for k in id_dict.keys()]),\
        ("One of these id keys is invalid: %s Valid keys are: %s."
         % (list(id_dict.keys()), valid_id_types))
    if 'doi' in id_dict.keys() and id_dict['doi'].lower().startswith('doi:'):
        id_dict['doi'] = id_dict['doi'][4:]
    content = None
    for id_type in valid_id_types:
        if id_type in id_dict.keys():
            content = download_article(id_dict[id_type], id_type)
            if content is not None:
                break
    else:
        logger.error("Could not download article with any of the ids: %s."
                     % str(id_dict))
    return content


def get_abstract(doi):
    """Get the abstract text of an article from Elsevier given a doi."""
    xml_string = download_article(doi)
    if xml_string is None:
        return None
    assert isinstance(xml_string, str)
    xml_tree = ET.XML(xml_string.encode('utf-8'), parser=UTB())
    if xml_tree is None:
        return None
    coredata = xml_tree.find('article:coredata', elsevier_ns)
    abstract = coredata.find('dc:description', elsevier_ns)
    abs_text = abstract.text
    return abs_text


def get_article(doi, output_format='txt'):
    """Get the full body of an article from Elsevier.

    Parameters
    ----------
    doi : str
        The doi for the desired article.
    output_format : 'txt' or 'xml'
        The desired format for the output. Selecting 'txt' (default) strips all
        xml tags and joins the pieces of text in the main text, while 'xml'
        simply takes the tag containing the body of the article and returns it
        as is . In the latter case, downstream code needs to be able to
        interpret Elsever's XML format.

    Returns
    -------
    content : str
        Either text content or xml, as described above, for the given doi.
    """
    xml_string = download_article(doi)
    if output_format == 'txt' and xml_string is not None:
        text = extract_text(xml_string)
        return text
    return xml_string


def extract_text(xml_string):
    """Get text from the body of the given Elsevier xml."""
    assert isinstance(xml_string, str)
    xml_tree = ET.XML(xml_string.encode('utf-8'), parser=UTB())
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
    return None


@lru_cache(maxsize=100)
@_ensure_api_keys('perform search')
def get_dois(query_str, count=100):
    """Search ScienceDirect through the API for articles.

    See http://api.elsevier.com/content/search/fields/scidir for constructing a
    query string to pass here.  Example: 'abstract(BRAF) AND all("colorectal
    cancer")'
    """
    url = '%s/%s' % (elsevier_search_url, query_str)
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


def get_piis(query_str):
    """Search ScienceDirect through the API for articles and return PIIs.

    Note that ScienceDirect has a limitation in which a maximum of 6,000
    PIIs can be retrieved for a given search and therefore this call is
    internally broken up into multiple queries by a range of years and the
    results are combined.

    Parameters
    ----------
    query_str : str
        The query string to search with

    Returns
    -------
    piis : list[str]
        The list of PIIs identifying the papers returned by the search
    """
    dates = range(1960, datetime.datetime.now().year)
    all_piis = flatten([get_piis_for_date(query_str, date) for date in dates])
    return all_piis


@lru_cache(maxsize=100)
@_ensure_api_keys('perform search')
def get_piis_for_date(query_str, date):
    """Search ScienceDirect with a query string constrained to a given year.

    Parameters
    ----------
    query_str : str
        The query string to search with
    date : str
        The year to constrain the search to

    Returns
    -------
    piis : list[str]
        The list of PIIs identifying the papers returned by the search
    """
    count = 200
    params = {'query': query_str,
              'count': count,
              'start': 0,
              'sort': '-coverdate',
              'date': date,
              'field': 'pii'}
    all_piis = []
    while True:
        res = requests.get(elsevier_search_url, params, headers=ELSEVIER_KEYS)
        if not res.status_code == 200:
            logger.info('Got status code: %d' % res.status_code)
            break
        res_json = res.json()
        entries = res_json['search-results']['entry']
        logger.info(res_json['search-results']['opensearch:totalResults'])
        if entries == [{'@_fa': 'true', 'error': 'Result set was empty'}]:
            logger.info('Search result was empty')
            return []
        piis = [entry['pii'] for entry in entries]
        all_piis += piis
        # Get next batch
        links = res_json['search-results'].get('link', [])
        cont = False
        for link in links:
            if link.get('@ref') == 'next':
                logger.info('Found link to next batch of results.')
                params['start'] += count
                cont = True
                break
        if not cont:
            break
    return all_piis


def download_from_search(query_str, folder):
    """Save raw text files based on a search for papers on ScienceDirect.

    This performs a search to get PIIs, downloads the XML corresponding to
    the PII, extracts the raw text and then saves the text into a file
    in the designated folder.

    Parameters
    ----------
    query_str : str
        The query string to search with
    folder : str
        The local path to an existing folder in which the text files
        will be dumped
    """
    piis = get_piis(query_str)
    for pii in piis:
        if os.path.exists(os.path.join(folder, '%s.txt' % pii)):
            continue
        logger.info('Downloading %s' % pii)
        xml = download_article(pii, 'pii')
        sleep(1)
        txt = extract_text(xml)
        if not txt:
            continue

        with open(os.path.join(folder, '%s.txt' % pii), 'wb') as fh:
            fh.write(txt.encode('utf-8'))


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
        return textwrap.fill(raw_text.text)
