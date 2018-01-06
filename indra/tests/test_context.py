from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import context_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr

@attr('webservice')
def test_get_protein_expression():
    res = context_client.get_protein_expression(['EGFR'], ['BT20_BREAST'])
    assert(res is not None)
    assert(res.get('BT20_BREAST') is not None)
    assert(res['BT20_BREAST'].get('EGFR') is not None)
    assert(res['BT20_BREAST']['EGFR'] > 1000)
    assert unicode_strs(res)


@attr('webservice')
def test_get_mutations():
    res = context_client.get_mutations(['BRAF'], ['A375_SKIN'])
    assert(res is not None)
    assert(res.get('A375_SKIN') is not None)
    assert(res['A375_SKIN'].get('BRAF') is not None)
    assert(res['A375_SKIN']['BRAF'] == ['V600E'])
    assert unicode_strs(res)


@attr('webservice')
def test_get_protein_expression_gene_missing():
    protein_amounts = context_client.get_protein_expression(['EGFR', 'XYZ'],
                                                            ['BT20_BREAST'])
    assert('BT20_BREAST' in protein_amounts)
    assert(protein_amounts['BT20_BREAST']['EGFR'] > 10000)
    assert(protein_amounts['BT20_BREAST']['XYZ'] is None)


@attr('webservice')
def test_get_protein_expression_cell_type_missing():
    protein_amounts = context_client.get_protein_expression(['EGFR'],
                                                            ['BT20_BREAST', 'XYZ'])
    assert('BT20_BREAST' in protein_amounts)
    assert(protein_amounts['BT20_BREAST']['EGFR'] > 10000)
    assert('XYZ' in protein_amounts)
    assert(protein_amounts['XYZ'] is None)


@attr('webservice')
def test_get_mutations_gene_missing():
    mutations = context_client.get_mutations(['BRAF', 'XYZ'], ['A375_SKIN'])
    assert('A375_SKIN' in mutations)
    assert(mutations['A375_SKIN']['BRAF'] == ['V600E'])
    assert(not mutations['A375_SKIN']['XYZ'])


@attr('webservice')
def test_get_mutations_cell_type_missing():
    mutations = context_client.get_mutations(['BRAF'], ['A375_SKIN', 'XYZ'])
    assert('A375_SKIN' in mutations)
    assert(mutations['A375_SKIN']['BRAF'] == ['V600E'])
    assert('XYZ' in mutations)
    assert(not mutations['XYZ']['BRAF'])
