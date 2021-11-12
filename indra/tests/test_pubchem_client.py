from indra.databases import pubchem_client


def test_get_inchi_key():
    ik = pubchem_client.get_inchi_key('5280613')
    assert ik == 'JLVSPVFPBBFMBE-HXSWCURESA-O', ik


def test_get_pmids():
    pmids = pubchem_client.get_pmids('2244')
    example_pmid = '19036898'
    assert example_pmid in pmids