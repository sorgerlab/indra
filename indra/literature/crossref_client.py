import requests
import json
from functools32 import lru_cache
import urllib
import re
from indra.literature import pubmed_client

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

def doi_query(pmid):
    # Get article metadata from PubMed
    article_xml = pubmed_client.get_article_xml(pmid)
    if article_xml is None:
        print "crossref_client: No article returned by PubMed."
        return None
    authors = [a.text for a in article_xml.findall('.//Author/LastName')]
    article_title = article_xml.find('ArticleTitle').text
    journal_title = article_xml.find('.//Journal/Title').text
    journal_issn = article_xml.find('.//Journal/ISSN').text
    if article_title is None or journal_title is None:
        print "crossref_client: Article or journal title missing."
        return None

    #combined_query = u'%s %s %s' % \
    #                 (article_title, journal_title, u' '.join(authors))
    combined_query = article_title

    url = crossref_search_url + 'dois?q=' + \
           urllib.quote_plus(combined_query.encode('UTF-8')) + \
          '&sort=score'
    res = requests.get(url)
    if res.status_code != 200:
        print 'Could not get DOI from CrossRef, code %d' % res.status_code
        return None
    raw_message = res.json()
    # Iterate over the search results, looking up XREF metadata
    for result_ix, result in enumerate(raw_message):
        if result_ix > 10:
            print "Giving up!"
            break
        doi_url = result['doi']
        metadata = get_metadata(doi_url)
        #journal_list = [r.upper()
        #                for r in metadata['container-title']]
        #xref_issn_list = metadata['ISSN']
        if metadata.get('author') is None:
            print "No author for XREF result, skipping"
            continue
        xref_author_list = [au['family'] for au in metadata['author']]
        author_matches = [pm_auth == xref_auth
                         for pm_auth, xref_auth in
                                zip(authors, xref_author_list)]
        if all(author_matches):
            print "ALl authors match, could be good!"
            m = re.match('^http://dx.doi.org/(.*)$', doi_url)
            doi = m.groups()[0]
            return doi
        elif len(xref_author_list) == len(authors):
            per_author_match = []
            for xref_auth, pm_auth in zip(xref_author_list, authors):
                if len(xref_auth) != len(pm_auth):
                    per_author_match.append(False)
                    continue
                char_matches = [c1 == c2 for c1, c2 in
                                zip(xref_auth, pm_auth)]
                if char_matches.count(False) <= 2:
                    per_author_match.append(True)
                else:
                    per_author_match.append(False)
            if all(per_author_match):
                print "Warning! Could be screwed up."
                print "xref", xref_author_list
                print "PM", authors
                m = re.match('^http://dx.doi.org/(.*)$', doi_url)
                doi = m.groups()[0]
                return doi
        else:
            print pmid, "authors don't match:"
            print "xref", xref_author_list
            print "PM", authors
        """
        if journal_issn not in xref_issn_list:
            print pmid, "XREF list", xref_issn_list
            print pmid, "pm journal ISSN", journal_issn
            print pmid, "Probably wrong"
        else:
            print pmid, "Maybe a match!"
            m = re.match('^http://dx.doi.org/(.*)$', doi_url)
            doi = m.groups()[0]
            return doi
        """
    return None

