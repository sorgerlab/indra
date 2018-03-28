from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import chebi_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr


def test_read_chebi_to_pubchem():
    (ctop, ptoc) = chebi_client._read_chebi_to_pubchem()
    assert ctop['85673'] == '91481662'
    assert ptoc['91481662'] == '85673'
    assert unicode_strs((ctop, ptoc))


def test_read_chebi_to_chembl():
    ctoc = chebi_client._read_chebi_to_chembl()
    assert ctoc['50729'] == 'CHEMBL58'
    assert unicode_strs(ctoc)


def test_cas_to_chebi():
    assert chebi_client.get_chebi_id_from_cas('23261-20-3') == '18035'
    assert chebi_client.get_chebi_id_from_cas('100-51-6') == '17987'
    assert chebi_client.get_chebi_id_from_cas('-1') is None
