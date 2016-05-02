from indra.literature import pubmed_client

def test_get_ids():
    ids = pubmed_client.get_ids('braf', retmax=10, db='pubmed')
    assert(len(ids) == 10)

def test_get_no_ids():
    ids = pubmed_client.get_ids('', retmax=10, db='pubmed')
    assert(not ids)

def test_get_pmc_ids():
    ids = pubmed_client.get_ids('braf', retmax=10, db='pmc')
    assert(len(ids) == 10)
    assert(len([i for i in ids if i.startswith('4') or
                i.startswith('3')]) == 10)

def test_get_abstract():
    abstract = pubmed_client.get_abstract('27085458')
    assert(abstract.startswith('Wilms'))
    assert(abstract.endswith('documented.'))

def test_get_no_abstract():
    abstract = pubmed_client.get_abstract('xx')
    assert(abstract is None)
