from indra.databases import doid_client
from indra.databases.doid_client import _client as client


def test_doid_client_loaded():
    assert 'doid' == client.prefix
    assert client.entries
    assert client.name_to_id


def test_doid_id_to_name():
    assert 'angiosarcoma' == \
        doid_client.get_doid_name_from_doid_id('DOID:0001816')


def test_doid_name_to_id():
    assert 'DOID:0001816' == \
        doid_client.get_doid_id_from_doid_name('angiosarcoma')


def test_doid_alt_id_to_doid_id():
    assert 'DOID:0001816' == \
        doid_client.get_doid_id_from_doid_alt_id('DOID:267')
