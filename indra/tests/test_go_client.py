from indra.databases import go_client


def test_invalid_id():
    go_name = go_client.get_go_label('34jkgfh')
    assert go_name is None


def test_go_id_lookup():
    go_id = 'GO:0001768'
    go_name = go_client.get_go_label(go_id)
    assert go_name == 'establishment of T cell polarity'


def test_go_label_to_id():
    assert go_client.get_go_id_from_label('mitochondrion inheritance') == \
        'GO:0000001'


def test_go_secondary_to_primary():
    assert go_client.get_primary_id('GO:0007067') == 'GO:0000278'


def test_get_valid_location():
    assert go_client.get_valid_location('0001669') == 'acrosomal vesicle'
    assert go_client.get_valid_location('GO:0001669') == 'acrosomal vesicle'
    assert go_client.get_valid_location('acrosomal vesicle') == \
        'acrosomal vesicle'
    assert go_client.get_valid_location('acrosome') == 'acrosomal vesicle'
