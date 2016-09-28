from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.literature import pubmed_client
from indra.util import unicode_strs

def test_get_ids():
    ids = pubmed_client.get_ids('braf', retmax=10, db='pubmed')
    assert(len(ids) == 10)
    assert unicode_strs(ids)

def test_get_no_ids():
    ids = pubmed_client.get_ids('', retmax=10, db='pubmed')
    assert(not ids)

def test_get_pmc_ids():
    ids = pubmed_client.get_ids('braf', retmax=10, db='pmc')
    assert(len(ids) == 10)
    assert(len([i for i in ids if i.startswith('5') or
                i.startswith('4')]) == 10)
    assert unicode_strs(ids)

def test_get_abstract():
    abstract = pubmed_client.get_abstract('27085458')
    assert(abstract.startswith('Wilms'))
    assert(abstract.endswith('documented.'))
    assert unicode_strs(abstract)
    # Try another one
    abstract = pubmed_client.get_abstract('27123883')
    assert unicode_strs(abstract)

def test_get_no_abstract():
    abstract = pubmed_client.get_abstract('xx')
    assert(abstract is None)

def test_get_ids_for_gene():
    ids = pubmed_client.get_ids_for_gene('EXOC1')
    assert ids
    assert unicode_strs(ids)

def test_get_metadata_for_ids():
    pmids = ['27123883', '27121204', '27115606']
    metadata = pubmed_client.get_metadata_for_ids(pmids)
    assert unicode_strs(metadata)
