from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import chebi_client
from indra.util import unicode_strs

def test_read_chebi_to_pubchem():
    (ctop, ptoc) = chebi_client.read_chebi_to_pubchem()
    assert ctop['85673'] == '252150010'
    assert ptoc['252150010'] == '85673'
    assert unicode_strs((ctop, ptoc))

def test_read_chebi_to_chembl():
    ctoc = chebi_client.read_chebi_to_chembl()
    assert ctoc['50729'] == 'CHEMBL58'
    assert unicode_strs(ctoc)


