from indra.databases import taxonomy_client


def test_name_lookup():
    assert taxonomy_client.get_taxonomy_id(
        'Severe acute respiratory syndrome coronavirus 2') == '2697049'
