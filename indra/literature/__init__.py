import urllib
from functools32 import lru_cache
from indra.literature import pubmed_client
from indra.literature import pmc_client
from indra.literature import crossref_client
from indra.literature import elsevier_client


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
