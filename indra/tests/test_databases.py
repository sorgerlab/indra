from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import get_identifiers_url


def test_chembl():
    cid = '1229517'
    assert get_identifiers_url('CHEMBL', cid) == \
        'https://identifiers.org/chembl.compound/CHEMBL%s' % cid
    assert get_identifiers_url('CHEMBL', 'CHEMBL%s' % cid) == \
        'https://identifiers.org/chembl.compound/CHEMBL%s' % cid


def test_signor():
    sid = 'SIGNOR-PF15'
    assert get_identifiers_url('SIGNOR', sid) == \
        'https://signor.uniroma2.it/relation_result.php?id=%s' % sid
