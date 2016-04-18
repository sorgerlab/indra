from indra.literature import crossref_client

test_doi = '10.1016/j.ccell.2016.02.010'

def test_get_metadata():
    metadata = crossref_client.get_metadata(test_doi)
    assert(metadata['DOI'] == test_doi)

    metadata = crossref_client.get_metadata('xyz')
    assert(metadata is None)

def test_get_publisher():
    publisher = crossref_client.get_publisher(test_doi)
    assert(publisher == 'Elsevier BV')

    publisher = crossref_client.get_publisher('xyz')
    assert(publisher is None)

def test_get_fulltext_links():
    links = crossref_client.get_fulltext_links(test_doi)
    assert(links[0]['content-type'] == 'text/xml')
    assert(links[1]['content-type'] == 'text/plain')

    links = crossref_client.get_fulltext_links('xyz')
    assert(links is None)

def test_get_license_links():
    links = crossref_client.get_license_links(test_doi)
    assert(links[0] == 'http://www.elsevier.com/tdm/userlicense/1.0/')

    links = crossref_client.get_license_links('xyz')
    assert(links is None)
