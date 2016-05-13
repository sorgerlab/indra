import urllib, urllib2
from functools32 import lru_cache
import xml.etree.ElementTree as ET
from indra.databases import hgnc_client

pubmed_search = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
pubmed_fetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
pmid_convert = 'http://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/'

@lru_cache(maxsize=100)
def send_request(url, data):
    try:
        req = urllib2.Request(url, data)
        res = urllib2.urlopen(req)
        xml_str = res.read()
        tree = ET.fromstring(xml_str)
    except:
        return None
    return tree


def get_ids(search_term, **kwargs):
    """Search Pubmed for paper IDs given a search term.

    The options are passed as named arguments. For details on parameters that
    can be used, see
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch Some useful
    parameters to pass are db='pmc' to search PMC instead of pubmed reldate=2
    to search for papers within the last 2 days mindate='2016/03/01',
    maxdate='2016/03/31' to search for papers in March 2016.
    """
    params = {'term': search_term,
              'retmax': 1000,
              'retstart': 0,
              'db': 'pubmed',
              'sort': 'pub+date'}
    for k, v in kwargs.iteritems():
        params[k] = v
    tree = send_request(pubmed_search, urllib.urlencode(params))
    if tree is None:
        return []
    if tree.find('ERROR') is not None:
        print tree.find('ERROR').text
        return []
    count = int(tree.find('Count').text)
    id_terms = tree.findall('IdList/Id')
    if id_terms is None:
        return []
    ids = [idt.text for idt in id_terms]
    if count != len(ids):
        print 'Not all ids were retrieved, limited at %d.' % params['retmax']
    return ids


def get_ids_for_gene(hgnc_name, **kwargs):
    """Get the curated set of articles for a gene in the Entrez database."""
    # Get the HGNC ID for the HGNC name
    hgnc_id = hgnc_client.get_hgnc_id(hgnc_name)
    if hgnc_id is None:
        raise ValueError('Invalid HGNC name.')
    # Get the Entrez ID
    entrez_id = hgnc_client.get_entrez_id(hgnc_id)
    if entrez_id is None:
        raise ValueError('Entrez ID not found in HGNC table.')
    # Query the Entrez Gene database
    params = {'db': 'gene',
              'retmode': 'xml',
              'id': entrez_id,
              }
    for k, v in kwargs.iteritems():
        params[k] = v
    tree = send_request(pubmed_fetch, urllib.urlencode(params))
    if tree is None:
        return []
    if tree.find('ERROR') is not None:
        print tree.find('ERROR').text
        return []
    # Get all PMIDs from the XML tree
    id_terms = tree.findall('.//PubMedId')
    if id_terms is None:
        return []
    # Use a set to remove duplicate IDs
    ids = list(set([idt.text for idt in id_terms]))
    return ids


def get_article_xml(pubmed_id):
    if pubmed_id.upper().startswith('PMID'):
        pubmed_id = pubmed_id[4:]
    params = {'db': 'pubmed',
                'retmode': 'xml',
                'id': pubmed_id}
    tree = send_request(pubmed_fetch, urllib.urlencode(params))
    if tree is None:
        return None
    article = tree.find('PubmedArticle/MedlineCitation/Article')
    return article # May be none


def get_title(pubmed_id):
    article = get_article_xml(pubmed_id)
    if article is None:
        return None
    title = article.find('ArticleTitle').text
    return title


def get_abstract(pubmed_id):
    article = get_article_xml(pubmed_id)
    if article is None:
        return None
    abstract = article.findall('Abstract/AbstractText')
    if abstract is None:
        return None
    else:
        abstract_text = ' '.join([' ' if abst.text is None 
                                    else abst.text for abst in abstract])
        return abstract_text

