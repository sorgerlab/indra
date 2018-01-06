from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import hgnc_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr


def test_get_uniprot_id():
    hgnc_id = '6840'
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    assert(uniprot_id == 'Q02750')
    assert unicode_strs(uniprot_id)


def test_get_uniprot_id_none():
    # This HGNC entry doesn't have a UniProt ID
    hgnc_id = '12027'
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    assert(uniprot_id is None)


def test_get_hgnc_name():
    hgnc_id = '3236'
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert(hgnc_name == 'EGFR')
    assert unicode_strs(hgnc_name)


@attr('webservice')
def test_get_hgnc_name_nonexistent():
    hgnc_id = '123456'
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert(hgnc_name is None)
    assert unicode_strs(hgnc_name)


def test_entrez_hgnc():
    entrez_id = '653509'
    hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
    assert(hgnc_id == '10798')


def test_entrez_hgnc_none():
    entrez_id = 'xxx'
    hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
    assert(hgnc_id is None)


def test_mouse_map():
    hgnc_id1 = hgnc_client.get_hgnc_from_mouse('109599')
    hgnc_id2 = hgnc_client.get_hgnc_from_mouse('MGI:109599')
    assert(hgnc_id1 == '4820')
    assert(hgnc_id2 == '4820')
    hgnc_id = hgnc_client.get_hgnc_from_mouse('xxx')
    assert(hgnc_id is None)


def test_rat_map():
    hgnc_id1 = hgnc_client.get_hgnc_from_rat('6496784')
    hgnc_id2 = hgnc_client.get_hgnc_from_rat('RGD:6496784')
    assert(hgnc_id1 == '44155')
    assert(hgnc_id2 == '44155')
    hgnc_id = hgnc_client.get_hgnc_from_rat('xxx')
    assert(hgnc_id is None)
