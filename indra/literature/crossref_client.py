import requests
import json
from functools32 import lru_cache

crossref_url = 'http://api.crossref.org/'

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
