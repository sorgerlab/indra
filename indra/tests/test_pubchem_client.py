from indra.databases import pubchem_client


def test_get_inchi_key():
    ik = pubchem_client.get_inchi_key('5280613')
    assert ik == 'JLVSPVFPBBFMBE-HXSWCURESA-O', ik


def test_get_pmids():
    pmids = pubchem_client.get_pmids('2244')
    example_pmid = '19036898'
    wrong_pmid = '17410596'
    assert example_pmid in pmids
    assert wrong_pmid not in pmids


def test_mesh_mappings():
    mesh_id = pubchem_client.get_mesh_id('56649450')  # Alpelisib
    assert mesh_id == 'C585539', mesh_id
    assert pubchem_client.get_mesh_id('abcd') is None


def test_get_substance_mesh_id():
    mesh_ids = pubchem_client.get_substance_mesh_id('27959613')
    example_mesh_id = 'D009570'
    wrong_mesh_id = 'D0074447'
    assert example_mesh_id in mesh_ids
    assert wrong_mesh_id not in mesh_ids
