from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import hgnc_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr


def test_get_uniprot_id():
    hgnc_id = '6840'
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    assert uniprot_id == 'Q02750'
    assert unicode_strs(uniprot_id)


def test_get_uniprot_id_none():
    # This HGNC entry doesn't have a UniProt ID
    hgnc_id = '37187'
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    assert uniprot_id is None, uniprot_id


def test_get_hgnc_name():
    hgnc_id = '3236'
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert hgnc_name == 'EGFR'
    assert unicode_strs(hgnc_name)


@attr('webservice')
def test_get_hgnc_name_nonexistent():
    hgnc_id = '123456'
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    assert hgnc_name is None
    assert unicode_strs(hgnc_name)


def test_entrez_hgnc():
    entrez_id = '653509'
    hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
    assert hgnc_id == '10798'


def test_entrez_hgnc_none():
    entrez_id = 'xxx'
    hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
    assert hgnc_id is None


def test_ensembl_hgnc():
    ensembl_id = 'ENSG00000006071'
    hgnc_id = hgnc_client.get_hgnc_from_ensembl(ensembl_id)
    assert hgnc_id == '59', hgnc_id
    assert hgnc_client.get_ensembl_id(hgnc_id) == ensembl_id


def test_mouse_map():
    hgnc_id1 = hgnc_client.get_hgnc_from_mouse('109599')
    hgnc_id2 = hgnc_client.get_hgnc_from_mouse('MGI:109599')
    assert hgnc_id1 == '4820'
    assert hgnc_id2 == '4820'
    hgnc_id = hgnc_client.get_hgnc_from_mouse('xxx')
    assert hgnc_id is None


def test_rat_map():
    hgnc_id1 = hgnc_client.get_hgnc_from_rat('6496784')
    hgnc_id2 = hgnc_client.get_hgnc_from_rat('RGD:6496784')
    assert hgnc_id1 == '44155'
    assert hgnc_id2 == '44155'
    hgnc_id = hgnc_client.get_hgnc_from_rat('xxx')
    assert hgnc_id is None


def test_is_category():
    assert hgnc_client.is_kinase('MAPK1')
    assert not hgnc_client.is_kinase('EGF')
    assert hgnc_client.is_phosphatase('PTEN')
    assert not hgnc_client.is_phosphatase('KRAS')
    assert hgnc_client.is_transcription_factor('FOXO3')
    assert not hgnc_client.is_transcription_factor('AKT1')


def test_get_current_id():
    # Current symbol
    assert hgnc_client.get_current_hgnc_id('BRAF') == '1097'
    # Outdated symbol, one ID
    assert hgnc_client.get_current_hgnc_id('SEPT7') == '1717'
    # Outdated symbol, multiple IDs
    ids = hgnc_client.get_current_hgnc_id('HOX1')
    assert len(ids) == 10
    assert '5101' in ids


def test_gene_type():
    assert hgnc_client.get_gene_type('1097') == 'gene with protein product'
    assert hgnc_client.get_gene_type('31547') == 'RNA, micro'
