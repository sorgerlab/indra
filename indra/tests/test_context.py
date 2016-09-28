from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import context_client
from indra.util import unicode_strs

def test_get_protein_expression():
    res = context_client.get_protein_expression('EGFR', 'BT20_BREAST')
    assert(res is not None)
    assert(res.get('EGFR') is not None)
    assert(res['EGFR'].get('BT20_BREAST') is not None)
    assert(res['EGFR']['BT20_BREAST'] > 1000)
    assert unicode_strs(res)

def test_get_mutations():
    res = context_client.get_mutations('BRAF', 'A375_SKIN')
    assert(res is not None)
    assert(res.get('BRAF') is not None)
    assert(res['BRAF'].get('A375_SKIN') is not None)
    assert(res['BRAF']['A375_SKIN'] == 1.0)
    assert unicode_strs(res)

