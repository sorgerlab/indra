from indra.databases import hp_client
from indra.databases.hp_client import _client as client


def test_hp_client_loaded():
    assert 'hp' == client.prefix
    assert client.entries
    assert client.name_to_id


def test_hp_id_to_name():
    assert 'Nocturia' == hp_client.get_hp_name_from_hp_id('HP:0000017')


def test_hp_name_to_id():
    assert 'HP:0000017' == hp_client.get_hp_id_from_hp_name('Nocturia')
