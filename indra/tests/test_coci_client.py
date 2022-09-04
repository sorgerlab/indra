from unittest import skip
from indra.literature import coci_client


@skip('COCI web service not working currently')
def test_citation_count():
    pmid = '24624335'
    doi = '10.1016/J.redox.2013.12.020'
    assert coci_client.get_citation_count_for_pmid(pmid) > 200
    assert coci_client.get_citation_count_for_doi(doi) > 200
