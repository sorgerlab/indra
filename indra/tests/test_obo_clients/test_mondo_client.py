from indra.databases import mondo_client

EXAMPLE_NAME = 'tenosynovial giant cell tumor, localized type'
EXAMPLE_ID = '0002399'
EXAMPLE_ALT_ID = '0024812'


def test_invalid_id():
    name = mondo_client.get_name_from_id('34jkgfh')
    assert name is None


def test_mondo_id_lookup():
    assert EXAMPLE_NAME == mondo_client.get_name_from_id(EXAMPLE_ID)


def test_mondo_label_to_id():
    assert EXAMPLE_ID == mondo_client.get_id_from_name(EXAMPLE_NAME)


def test_mondo_secondary_to_primary():
    assert EXAMPLE_ID == mondo_client.get_id_from_alt_id(EXAMPLE_ALT_ID)
