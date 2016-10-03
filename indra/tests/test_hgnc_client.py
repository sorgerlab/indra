from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import hgnc_client
from indra.util import unicode_strs

def test_get_uniprot_id():
    hgnc_id = '6840'
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    assert(uniprot_id == 'Q02750')
    assert unicode_strs(uniprot_id)

def test_get_hgnc_name():
    hgnc_id = '3236'
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert(hgnc_name == 'EGFR')
    assert unicode_strs(hgnc_name)

def test_get_hgnc_name_nonexistent():
    hgnc_id = '123456'
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert(hgnc_name is None)
    assert unicode_strs(hgnc_name)

def test_entrez_hgnc():
    entrez_id = '653509'
    hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
    assert(hgnc_id == '10798')
