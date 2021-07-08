from indra.databases import mondo_client

EXAMPLE_NAME = 'tenosynovial giant cell tumor, localized type'
EXAMPLE_ID = '0002399'
EXAMPLE_ALT_ID = '0024812'


def test_invalid_id():
    name = mondo_client.get_name_from_id('34jkgfh')
    assert name is None


def test_mondo_id_lookup():
    name = mondo_client.get_name_from_id(EXAMPLE_ID)
    assert name is not None
    assert EXAMPLE_NAME == name


def test_mondo_label_to_id():
    identifier = mondo_client.get_id_from_name(EXAMPLE_NAME)
    assert identifier is not None
    assert EXAMPLE_ID == identifier


def test_mondo_secondary_to_primary():
    identifier = mondo_client.get_id_from_alt_id('MONDO:0018220')
    assert identifier is not None
    assert '0002413' == identifier
