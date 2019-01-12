import logging
from collections import defaultdict

from indra.literature.elsevier_client import logger as elsevier_logger
from indra.literature import pubmed_client, pmc_client, elsevier_client

logger = logging.getLogger('disambiguate')

# the elsevier_client will log messages that it is safe to ignore
elsevier_logger.setLevel(logging.WARNING)


def get_fulltexts_from_entrez(hgnc_name):
    pmids = pubmed_client.get_ids_for_gene(hgnc_name)
    articles = (pubmed_client.get_article_xml(pmid) for pmid in pmids)
    fulltexts = [_universal_extract_text(article) for article in articles]
    return fulltexts


def _universal_extract_text(xml):
    # first try to parse the xml as if it came from elsevier. if we do not
    # have valid elsevier xml this will throw an exception.
    # the text extraction function in the pmc client may not throw an
    # exception when parsing elsevier xml, silently processing the xml
    # incorrectly
    try:
        fulltext = elsevier_client.extract_text(xml)
    except Exception:
        try:
            fulltext = pmc_client.extract_text(xml)
        except Exception:
            # fall back by returning input string unmodified
            fulltext = xml
    return fulltext


def _get_text_from_pmids(pmids):
    pmc_content = set(pubmed_client.filter_pmids(pmids))
    pmc_ids = (pmc_client.id_lookup(pmid, idtype='pmid')['pmcid']
               for pmid in pmc_content)
    pmc_xmls = (pmc_client.get_xml(pmc_id) for pmc_id in pmc_ids)
    pmc_texts = set(_universal_extract_text(xml) for xml in pmc_xmls)

    other_content = set(pmids) - pmc_content
    ids = (pmc_client.id_lookup(pmid, idtype='pmid') for pmid in pmids)
    elsevier_content = (elsevier_client.download_article_from_id(pmid)
                        for pmid in pmids)
    
