from indra.databases import biolookup_client


def test_lookup_curie():
    curie = 'pubchem.compound:40976'
    res = biolookup_client.lookup_curie(curie)
    assert res['name'] == '(17alpha)-13-Ethyl-17-hydroxy-11-methylene' \
        '-18,19-dinorpregn-4-en-20-yn-3-one', res


def test_lookup():
    res = biolookup_client.lookup('FPLX', 'ERK')
    assert res['name'] == 'ERK', res


def test_get_name():
    res = biolookup_client.get_name('CHEBI', 'CHEBI:408174')
    assert res == 'arformoterol', res