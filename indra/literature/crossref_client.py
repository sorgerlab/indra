import requests
import json
from functools32 import lru_cache
import urllib
import re

crossref_url = 'http://api.crossref.org/'
crossref_search_url = 'http://search.crossref.org/'

@lru_cache(maxsize=100)
def get_metadata(doi):
    """Returns the metadata of an article given its DOI from CrossRef
    as a JSON dict"""
    url = crossref_url + 'works/' + doi
    res = requests.get(url)
    if res.status_code != 200:
        print 'Could not get CrossRef metadata, code %d' % res.status_code
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

def get_license_links(doi):
    metadata = get_metadata(doi)
    if metadata is None:
        return None
    licenses = metadata.get('license')
    if licenses is None:
        return None
    urls = [l.get('URL') for l in licenses]
    return urls

@lru_cache(maxsize=100)
def doi_query(title):
    # If None or empty string return None
    if not title:
        return None
    url = crossref_search_url + 'dois?q=' + \
          urllib.quote_plus(title.encode('UTF-8')) + \
          'sort=score'
    res = requests.get(url)
    if res.status_code != 200:
        print 'Could not get DOI from CrossRef, code %d' % res.status_code
        return None
    raw_message = res.json()
    first_result = raw_message[0]
    doi_url = first_result['doi'] # Return DOI of first result (Yikes!)
    m = re.match('^http://dx.doi.org/(.*)$', doi_url)
    doi = m.groups()[0]
    return doi
