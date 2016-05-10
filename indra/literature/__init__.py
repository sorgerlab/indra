import urllib
from functools32 import lru_cache
from indra.literature import pubmed_client
from indra.literature import pmc_client
from indra.literature import crossref_client
from indra.literature import elsevier_client


def id_lookup(paper_id, idtype):
    """Take an ID of type PMID, PMCID, or DOI and lookup the other IDs.

    If the DOI is not found in Pubmed, try to obtain the DOI by doing a
    reverse-lookup of the DOI in CrossRef using article metadata.

    Parameters
    ----------
    paper_id : string
        ID of the article.
    idtype : 'pmid', 'pmcid', or 'doi
        Type of the ID.
    """
    if idtype not in ('pmid', 'pmcid', 'doi'):
        raise ValueError("Invalid idtype %s; must be 'pmid', 'pmcid', "
                         "or 'doi'." % idtype)

    ids = {'doi': None, 'pmid': None, 'pmcid': None}
    pmc_id_results = pmc_client.id_lookup(paper_id, idtype)
    ids['pmid'] = pmc_id_results.get('pmid')
    ids['pmcid'] = pmc_id_results.get('pmcid')
    ids['doi'] = pmc_id_results.get('doi')
    # If we gave a DOI, then our work is done after looking for PMID and PMCID
    if idtype == 'doi':
        return ids
    # If we gave a PMID or PMCID, we need to check to see if we got a DOI.
    # If we got a DOI back, we're done.
    elif ids.get('doi'):
        return ids
    # If we get here, then we've given PMID or PMCID and don't have a DOI yet.
    # If we gave a PMCID and have neither a PMID nor a DOI, then we'll run
    # into problems later on when we try to the reverse lookup using CrossRef.
    # So we bail here and return what we have (PMCID only) with a warning.
    if ids.get('pmcid') and ids.get('doi') is None and ids.get('pmid') is None:
        warnings.warn('%s: PMCID without PMID or DOI' % ids.get('pmcid'))
        return ids
    # To clarify the state of things at this point:
    assert ids.get('pmid') is not None
    assert ids.get('doi') is None
    # Now we try to get the DOI from CrossRef:
    ids['doi'] = crossref_client.doi_query(pmid)
    # It may still be None, but at this point there's nothing we can do...
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
