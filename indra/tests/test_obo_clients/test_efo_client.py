from indra.databases import efo_client
from indra.databases.efo_client import _client as client


def test_efo_client_loaded():
    assert 'efo' == client.prefix
    assert client.entries
    assert client.name_to_id


def test_efo_id_to_name():
    assert 'muscle measurement' == \
        efo_client.get_efo_name_from_efo_id('0004515')


def test_efo_name_to_id():
    assert '0004515' == \
        efo_client.get_efo_id_from_efo_name('muscle measurement')
