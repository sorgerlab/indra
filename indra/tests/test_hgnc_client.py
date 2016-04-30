from indra.databases import hgnc_client

def test_get_uniprot_id():
    hgnc_id = '6840'
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    assert(uniprot_id == 'Q02750')

def test_get_hgnc_name():
    hgnc_id = '3236'
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert(hgnc_name == 'EGFR')

def test_get_hgnc_name_nonexistent():
    hgnc_id = '123456'
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert(hgnc_name is None)
