import urllib, urllib2
from functools32 import lru_cache
import xml.etree.ElementTree as ET

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

def get_ids(search_term, retmax=1000, db='pubmed'):
    params = {'db': db,
                'term': search_term,
                'sort': 'pub+date',
                'retstart': 0,
                'retmax': retmax}
    tree = send_request(pubmed_search, urllib.urlencode(params))
    if tree is None:
        return []
    count = int(tree.find('Count').text)
    id_terms = tree.findall('IdList/Id')
    if id_terms is None:
        return []
    ids = [idt.text for idt in id_terms]
    if count != len(ids):
        print 'Not all ids were retrieved, limited at %d.' % retmax
    return ids

def get_abstract(pubmed_id):
    if pubmed_id.upper().startswith('PMID'):
        pubmed_id = pubmed_id[4:]
    params = {'db': 'pubmed',
                'retmode': 'xml',
                'rettype': 'abstract',
                'id': pubmed_id}
    tree = send_request(pubmed_fetch, urllib.urlencode(params))
    if tree is None:
        return None
    article = tree.find('PubmedArticle/MedlineCitation/Article')
    if article is None:
        return None
    abstract = article.findall('Abstract/AbstractText')
    if abstract is None:
        return None
    else:
        abstract_text = ' '.join([' ' if abst.text is None 
                                    else abst.text for abst in abstract])
        return abstract_text

def pmid_to_doi(pubmed_id):
    if pubmed_id.upper().startswith('PMID'):
        pubmed_id = pubmed_id[4:]
    url = pmid_convert
    data = {'ids': pubmed_id}
    tree = send_request(url, urllib.urlencode(data))
    if tree is None:
        return None
    record = tree.find('record')
    if record is None:
        return None
    doi = record.attrib.get('doi')
    return doi 
