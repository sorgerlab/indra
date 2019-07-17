from indra.databases import mirbase_client


def test_mirbase_id_to_name():
    assert 'hsa-mir-19b-2' == \
            mirbase_client.get_mirbase_name_from_mirbase_id('MI0000075')


def test_mirbase_name_to_id():
    assert 'MI0000075' == \
            mirbase_client.get_mirbase_id_from_mirbase_name('hsa-mir-19b-2')


def test_mirbase_hgnc_mappings():
    assert '31476' == mirbase_client.get_hgnc_id_from_mirbase_id('MI0000060')
    assert mirbase_client.get_hgnc_id_from_mirbase_id('MI0000056') is None, \
                                            'This one is from C. elegans'


def test_hgnc_mirbase_mappings():
    assert 'MI0000060' == mirbase_client.get_mirbase_id_from_hgnc_id('31476')
    assert mirbase_client.get_mirbase_id_from_hgnc_id('6893') is None, \
                                            'HGNC:6893!MAPT is a protein'


def test_hgnc_symbol_mirbase_mappings():
    assert 'MI0000075' == \
                    mirbase_client.get_mirbase_id_from_hgnc_symbol('MIR19B2')
    assert mirbase_client.get_mirbase_id_from_hgnc_symbol('MAPT') is None, \
                                            'HGNC:6893!MAPT is a protein'
