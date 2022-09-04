from indra.databases import biolookup_client


def test_lookup_curie():
    curie = 'pubchem.compound:40976'
    res = biolookup_client.lookup_curie(curie)
    assert res['name'] == '(17R)-13-ethyl-17-ethynyl-17-hydroxy-11-' \
        'methylidene-2,6,7,8,9,10,12,14,15,16-decahydro-1H-' \
        'cyclopenta[a]phenanthren-3-one', res


def test_lookup():
    res = biolookup_client.lookup('HGNC', '1097')
    assert res['name'] == 'BRAF', res


def test_get_name():
    res = biolookup_client.get_name('CHEBI', 'CHEBI:408174')
    assert res == 'arformoterol', res
