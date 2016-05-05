import urllib
from functools32 import lru_cache
from indra.literature import pubmed_client
from indra.literature import pmc_client
from indra.literature import crossref_client
from indra.literature import elsevier_client

@lru_cache(maxsize=100)
def id_lookup(paper_id):
    """This function takes an ID of the form
    PMID*, PMC* or DOI*, uses the Pubmed ID mapping
    service and looks up all other IDs from one
    of these. The IDs are returned in a dictionary"""
    if paper_id.upper().startswith('PMID'):
        idtype = 'pmid'
        paper_id = paper_id[4:]
    elif paper_id.upper().startswith('PMC'):
        idtype = 'pmcid'
        paper_id = paper_id[3:]
    elif paper_id.upper().startswith('DOI'):
        idtype = 'doi'
        paper_id = paper_id[3:]
    else:
        idtype = None
    url = pubmed_client.pmid_convert
    data = {'ids': paper_id}
    if idtype is not None:
        data['idtype'] = idtype
    tree = pubmed_client.send_request(url, urllib.urlencode(data))
    if tree is None:
        return {}
    record = tree.find('record')
    if record is None:
        return {}
    doi = record.attrib.get('doi')
    pmid = record.attrib.get('pmid')
    pmcid = record.attrib.get('pmcid')
    # If we've failed to find the DOI, and we have the PMID, we do a
    # fallback--get the title from PubMed and then search CrossRef for the DOI
    if doi is None and pmid is not None:
        print "Couldn't find the DOI in PubMed, querying CrossRef"
        title = pubmed_client.get_title(pmid)
        doi = crossref_client.doi_query(title)
        print "Found doi: %s" % doi
    ids = {'doi': doi,
           'pmid': pmid,
           'pmcid': pmcid}
    return ids

def get_full_text(paper_id):
    ids = id_lookup(paper_id)
    pmcid = ids.get('pmcid')
    pmid = ids.get('pmid')
    doi = ids.get('doi')
    # First try to find paper via PMC
    if pmcid:
        nxml = pmc_client.get_xml(pmcid)
        if nxml:
            return nxml, 'nxml'

    # If it does not have PMC NXML then we need to use
    # a dedicated literature client
    if doi:
        publisher = crossref_client.get_publisher(doi)
        links = crossref_client.get_fulltext_links(doi)
        if publisher == 'Elsevier BV':
            full_text = elsevier_client.get_article(doi, output='txt')
            if full_text:
                return full_text, 'txt'

    # If no full text resource is available then return the 
    # abstract
    if pmid:
        abstract = pubmed_client.get_abstract(pmid)
        return abstract, 'abstract'

    return None, None
