"""
For information on the Elsevier API, see:
  - API Specification: http://dev.elsevier.com/api_docs.html
  - Authentication: https://dev.elsevier.com/tecdoc_api_authentication.html
"""
import os
import re
import logging
import textwrap
import datetime
import xml.etree.ElementTree as ET
import requests
from time import sleep
from indra.util import flatten
from indra import has_config, get_config
from functools import lru_cache, wraps
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger(__name__)


# THE ELSEVIER API URL: ***MUST BE HTTPS FOR SECURITY***
elsevier_api_url = 'https://api.elsevier.com/content'  # <--- HTTPS
elsevier_article_url_fmt = '%s/article/%%s' % elsevier_api_url
elsevier_search_url = '%s/search/sciencedirect' % elsevier_api_url
elsevier_entitlement_url = '%s/article/entitlement/doi' % elsevier_api_url

# Namespaces for Elsevier XML elements
elsevier_ns = {'dc': 'http://purl.org/dc/elements/1.1/',
               'article': 'http://www.elsevier.com/xml/svapi/article/dtd',
               'ja': 'http://www.elsevier.com/xml/ja/dtd',
               'xocs': 'http://www.elsevier.com/xml/xocs/dtd',
               'common': 'http://www.elsevier.com/xml/common/dtd',
               'atom': 'http://www.w3.org/2005/Atom',
               'prism': 'http://prismstandard.org/namespaces/basic/2.0/',
               'book': 'http://www.elsevier.com/xml/bk/dtd'}
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
        You can enter any combination of eid, doi, pmid, and/or pii. Ids will
        be checked in that order, until either content has been found or all
        ids have been checked.

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
    paragraphs = extract_paragraphs(xml_string)
    if paragraphs:
        return '\n'.join(re.sub('\s+', ' ', p) for p in paragraphs) + '\n'
    else:
        return None


def extract_paragraphs(xml_string):
    """Get paragraphs from the body of the given Elsevier xml."""
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
        return [raw_text]
    return None


@lru_cache(maxsize=100)
@_ensure_api_keys('perform search')
def get_dois(query_str, year=None, loaded_after=None):
    """Search ScienceDirect through the API for articles and return DOIs.

    Parameters
    ----------
    query_str : str
        The query string to search with.
    year : Optional[str]
        The year to constrain the search to.
    loaded_after : Optional[str]
        Date formatted as 'yyyy-MM-dd'T'HH:mm:ssX' to constrain the search
        to articles loaded after this date. Example: 2019-06-01T00:00:00Z

    Returns
    -------
    dois : list[str]
        The list of DOIs identifying the papers returned by the search.
    """
    dois = search_science_direct(
        query_str, field_name='doi', year=year, loaded_after=loaded_after)
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
def get_piis_for_date(query_str, year=None, loaded_after=None):
    """Search ScienceDirect through the API for articles and return PIIs.

    Parameters
    ----------
    query_str : str
        The query string to search with.
    year : Optional[str]
        The year to constrain the search to.
    loaded_after : Optional[str]
        Date formatted as 'yyyy-MM-dd'T'HH:mm:ssX' to constrain the search
        to articles loaded after this date. Example: 2019-06-01T00:00:00Z

    Returns
    -------
    piis : list[str]
        The list of PIIs identifying the papers returned by the search.
    """
    piis = search_science_direct(
        query_str, field_name='pii', year=year, loaded_after=loaded_after)
    return piis


@lru_cache(maxsize=100)
@_ensure_api_keys('perform search')
def search_science_direct(query_str, field_name, year=None, loaded_after=None):
    """Search ScienceDirect for a given field with a query string.

    Users can specify which field they are interested in and only values from
    that field will be returned. It is also possible to restrict the search
    either to a specific year of publication or to papers published after a
    specific date.

    Parameters
    ----------
    query_str : str
        The query string to search with.
    field_name : str
        A name of the field of interest to be returned. Accepted values are:
        authors, doi, loadDate, openAccess, pages, pii, publicationDate,
        sourceTitle, title, uri, volumeIssue.
    year : Optional[str]
        The year to constrain the search to.
    loaded_after : Optional[str]
        Date formatted as 'yyyy-MM-dd'T'HH:mm:ssX' to constrain the search
        to articles loaded after this date.

    Returns
    -------
    all_parts : list[str]
        The list of values from the field of interest identifying the papers
        returned by the search.
    """
    count = 100
    params = {'qs': query_str,
              'display': {
                  'offset': 0,
                  'show': count,
                  'sortBy': 'date'},
              'field': 'pii'}
    if year:
        params['date'] = year
    if loaded_after:
        params['loadedAfter'] = loaded_after
    all_parts = []
    while True:
        res = requests.put(
            elsevier_search_url, json=params, headers=ELSEVIER_KEYS)
        if not res.status_code == 200:
            logger.info('Got status code: %d' % res.status_code)
            break
        res_json = res.json()
        total_results = res_json['resultsFound']
        if total_results == 0:
            logger.info('Search result was empty')
            return []
        entries = res_json['results']
        parts = [entry[field_name] for entry in entries]
        all_parts += parts
        # Get next batch
        cont = False
        # We can only set offset up to 6000
        if (params['display']['offset'] + count) <= min(total_results, 6000):
            params['display']['offset'] += count
            cont = True
            # There is a quota on number of requests, wait to continue
            sleep(1)
        if not cont:
            break
    return all_parts


def download_from_search(query_str, folder, do_extract_text=True,
                         max_results=None):
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
    do_extract_text : bool
        Choose whether to extract text from the xml, or simply save the raw xml
        files. Default is True, so text is extracted.
    max_results : int or None
        Default is None. If specified, limit the number of results to the given
        maximum.
    """
    piis = get_piis(query_str)
    for pii in piis[:max_results]:
        if os.path.exists(os.path.join(folder, '%s.txt' % pii)):
            continue
        logger.info('Downloading %s' % pii)
        xml = download_article(pii, 'pii')
        sleep(1)
        if do_extract_text:
            txt = extract_text(xml)
            if not txt:
                continue

            with open(os.path.join(folder, '%s.txt' % pii), 'wb') as fh:
                fh.write(txt.encode('utf-8'))
        else:
            with open(os.path.join(folder, '%s.xml' % pii), 'wb') as fh:
                fh.write(xml.encode('utf-8'))
    return


def _get_article_body(full_text_elem):
    possible_paths = [
        'xocs:doc/xocs:serial-item/ja:article/ja:body',
        'xocs:doc/xocs:serial-item/ja:simple-article/ja:body',
        'xocs:doc/xocs:serial-item/ja:converted-article/ja:body',
        'xocs:doc/xocs:nonserial-item/book:chapter',
        'xocs:doc/xocs:nonserial-item/book:fb-non-chapter'
        ]
    for pth in possible_paths:
        main_body = full_text_elem.find(pth, elsevier_ns)
        if main_body is not None:
            logger.info("Found main body element: \"%s\"" % pth)
            return _get_sections(main_body)
        logger.info("Could not find main body element: \"%s\"." % pth)
    return None


def _get_sections(main_body_elem):
    # Get content sections
    possible_paths = ['common:sections/common:section', 'common:section',
                      'common:sections']
    for pth in possible_paths:
        sections = main_body_elem.findall(pth, elsevier_ns)
        if len(sections):
            logger.info("Found sections in main body using \"%s\"" % pth)
            break
        logger.info("Found no sections in main body with \"%s\"" % pth)
    else:
        return None

    # Concatenate the section content
    paragraphs = []
    for s in sections:
        # Paragraphs that are directly under the section
        pars = s.findall('common:para', elsevier_ns)
        # Paragraphs that are under a section within the section
        pars += s.findall('common:section/common:para', elsevier_ns)
        for p in pars:
            content = ' '.join(p.itertext())
            paragraphs.append(content)
    return paragraphs


def _get_raw_text(full_text_elem):
    # Look for raw_text
    raw_text = full_text_elem.find('xocs:doc/xocs:rawtext', elsevier_ns)
    if raw_text is None:
        logger.info("Could not find rawtext element xocs:doc/xocs:rawtext")
        return None
    else:
        logger.info("Found rawtext element xocs:doc/xocs:rawtext")
        return textwrap.fill(raw_text.text)
