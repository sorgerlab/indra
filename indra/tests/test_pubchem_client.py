from indra.databases import pubchem_client


def test_get_inchi_key():
    ik = pubchem_client.get_inchi_key('5280613')
    assert ik == 'JLVSPVFPBBFMBE-HXSWCURESA-O', ik
