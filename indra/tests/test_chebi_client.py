from indra.databases import chebi_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr


def test_read_chebi_to_pubchem():
    (ctop, ptoc) = chebi_client._read_chebi_to_pubchem()
    assert ctop['CHEBI:85673'] == '11652416'
    assert ptoc['11652416'] == 'CHEBI:85673'
    assert unicode_strs((ctop, ptoc))


def test_chebi_pubchem_mapping():
    # This is a non-trivial mapping since there are multiple mappings
    # reported by ChEBI and we need to choose the right one based on
    # InChIKey matches.
    assert chebi_client.get_chebi_id_from_pubchem('5287993') == 'CHEBI:3528'
    assert chebi_client.get_pubchem_id('CHEBI:3528') == '5287993'


def test_read_chebi_to_chembl():
    chebi_to_chembl, _ = chebi_client._read_chebi_to_chembl()
    assert chebi_to_chembl['CHEBI:50729'] == 'CHEMBL58'
    assert unicode_strs(chebi_to_chembl)


def test_chebi_chembl():
    assert chebi_client.get_chebi_id_from_chembl('CHEMBL525191') == \
        'CHEBI:83405'
    assert chebi_client.get_chembl_id('CHEBI:83405') == 'CHEMBL525191'


def test_cas_to_chebi():
    assert chebi_client.get_chebi_id_from_cas('23261-20-3') == 'CHEBI:18035'
    assert chebi_client.get_chebi_id_from_cas('100-51-6') == 'CHEBI:17987'
    assert chebi_client.get_chebi_id_from_cas('-1') is None


def test_chebi_id_to_name():
    name = chebi_client.get_chebi_name_from_id('CHEBI:63637')
    assert name == 'vemurafenib', name


def test_chebi_name_to_id():
    cid = chebi_client.get_chebi_id_from_name('vemurafenib')
    assert cid == 'CHEBI:63637', cid


@attr('webservice')
def test_chebi_name_from_web():
    name = chebi_client.get_chebi_name_from_id_web('63637')
    assert name == 'vemurafenib'
    name = chebi_client.get_chebi_name_from_id_web('44215')
    assert name == 'NAD zwitterion'


@attr('webservice')
def test_inchi_key():
    ik = chebi_client.get_inchi_key('2150')
    assert ik == 'NVKAWKQGWWIWPM-MISPCMORSA-N'


def test_hmdb_to_chebi():
    chebi_id = chebi_client.get_chebi_id_from_hmdb('HMDB0000122')
    assert chebi_id == 'CHEBI:15903', chebi_id


def test_chebi_to_primary():
    assert chebi_client.get_primary_id('CHEBI:6281') == 'CHEBI:17490'
    assert chebi_client.get_primary_id('CHEBI:161680') == 'CHEBI:161680'
