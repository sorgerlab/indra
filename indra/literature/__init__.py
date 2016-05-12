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
    ids['doi'] = crossref_client.doi_query(ids['pmid'])
    # It may still be None, but at this point there's nothing we can do...
    return ids


def get_full_text(paper_id, idtype):
    ids = id_lookup(paper_id, idtype)
    pmcid = ids.get('pmcid')
    pmid = ids.get('pmid')
    doi = ids.get('doi')
    # First try to find paper via PMC
    if pmcid:
        nxml = pmc_client.get_xml(pmcid)
        if nxml:
            return nxml, 'nxml'
    # If we got here, it means we didn't find the full text in PMC, so we'll
    # need either the DOI (for lookup in CrossRef) and/or the PMID (so we
    # can fall back on the abstract. If by some strange turn we have neither,
    # give up now.
    if not doi and not pmid:
        return None, None

    # If it does not have PMC NXML then we attempt to obtain the full-text
    # through the CrossRef Click-through API
    if doi:
        # Check if there are any full text links
        links = crossref_client.get_fulltext_links(doi)
        # Get publisher
        publisher = crossref_client.get_publisher(doi)
        if links:
            # Utility function to get particular links by content-type
            def lookup_content_type(link_list, content_type):
                content_list = [l.get('URL') for l in link_list
                                if l.get('content-type') == content_type]
                return None if not content_list else content_list[0]
            # Check for XML first
            if lookup_content_type(links, 'text/xml'):
                pass
            elif lookup_content_type(links, 'text/plain'):
                pass
            elif lookup_content_type(links, 'application/pdf'):
                pass
            elif lookup_content_type(links, 'unknown'):
                pass
            else:
                raise Exception("Unknown content type(s): %s" % links)
        # No full text links :( Check and see if the publisher is one we have
        # a web scraper for
        elif publisher in []:
            pass
        # No full text links and not a publisher we support. We'll have to
        # fall back to the abstract.
        elif pmid:
            abstract = pubmed_client.get_abstract(pmid)
            return abstract, 'abstract'
        # We have a useless DOI and no PMID. Give up.
        else:
            return None, None
    # We don't have a DOI but we're guaranteed to have a PMID at this point,
    # so we fall back to the abstract:
    else:
        abstract = pubmed_client.get_abstract(pmid)
        return abstract, 'abstract'

    # We'll only get here if we've missed a combination of conditions
    assert False

