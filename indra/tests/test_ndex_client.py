from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import ndex_client


def test_ndex_ver():
    assert ndex_client._increment_ndex_ver(None) == '1.0'
    assert ndex_client._increment_ndex_ver('') == '1.0'
    assert ndex_client._increment_ndex_ver('1.0') == '1.1'
    assert ndex_client._increment_ndex_ver('1.9') == '1.10'
    assert ndex_client._increment_ndex_ver('2.10') == '2.11'

