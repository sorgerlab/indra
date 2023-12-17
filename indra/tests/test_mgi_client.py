from indra.databases import mgi_client


def test_lookups():
    assert mgi_client.get_id_from_name('Braf') == '88190'
    assert mgi_client.get_name_from_id('1926283') == 'Pgap6'
    assert mgi_client.get_id_from_name_synonym('Pgap6') == '1926283'
    assert mgi_client.get_id_from_name_synonym('Tmem8') == '1926283'
    assert isinstance(mgi_client.get_id_from_name_synonym('EGF-TM7'), list)
    assert mgi_client.get_ensembl_id('88190') == 'ENSMUSG00000002413'
