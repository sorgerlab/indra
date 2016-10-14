from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.literature import crossref_client
from indra.util import unicode_strs

test_doi = '10.1016/j.ccell.2016.02.010'

example_ids = {'pmid': '25361007',
               'pmcid': 'PMC4322985',
               'doi': '10.18632/oncotarget.2555'}
def test_doi_query():
    mapped_doi = crossref_client.doi_query(example_ids['pmid'])
    assert mapped_doi == example_ids['doi']
    assert unicode_strs(mapped_doi)

def test_get_metadata():
    metadata = crossref_client.get_metadata(test_doi)
    assert(metadata['DOI'] == test_doi)
    assert unicode_strs(metadata)
    metadata = crossref_client.get_metadata('xyz')
    assert(metadata is None)

def test_get_publisher():
    publisher = crossref_client.get_publisher(test_doi)
    assert(publisher == 'Elsevier BV')
    assert unicode_strs(publisher)
    publisher = crossref_client.get_publisher('xyz')
    assert(publisher is None)

def test_get_fulltext_links():
    links = crossref_client.get_fulltext_links(test_doi)
    assert(links[0]['content-type'] == 'text/xml')
    assert(links[1]['content-type'] == 'text/plain')
    assert unicode_strs(links)
    links = crossref_client.get_fulltext_links('xyz')
    assert(links is None)

def test_get_license_links():
    links = crossref_client.get_license_links(test_doi)
    assert(links[0] == 'http://www.elsevier.com/tdm/userlicense/1.0/')
    assert unicode_strs(links)
    links = crossref_client.get_license_links('xyz')
    assert(links is None)

def test_get_url():
    url = crossref_client.get_url(test_doi)
    assert url == 'http://dx.doi.org/10.1016/j.ccell.2016.02.010'
    assert unicode_strs(url)
    url = crossref_client.get_url('xyz')
    assert url is None
