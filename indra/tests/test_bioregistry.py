from indra.databases import bioregistry


def test_get_ns_from_bioregistry():
    assert bioregistry.get_ns_from_bioregistry('xxxx') is None
    assert bioregistry.get_ns_from_bioregistry('noncodev4.rna') == 'NONCODE'
    assert bioregistry.get_ns_from_bioregistry('chebi') == 'CHEBI'


def test_get_ns_id_from_bioregistry():
    assert bioregistry.get_ns_id_from_bioregistry('xxxx', 'xxxx') == \
        (None, None)
    assert bioregistry.get_ns_id_from_bioregistry('chebi', '3696') == \
        ('CHEBI', 'CHEBI:3696')
    assert bioregistry.get_ns_id_from_bioregistry('hgnc', '1097') == \
        ('HGNC', '1097')


def test_get_ns_id_from_bioregistry_curie():
    assert bioregistry.get_ns_id_from_bioregistry_curie('xxxx:xxxx') == \
        (None, None)
    assert bioregistry.get_ns_id_from_bioregistry_curie('chebi:3696') == \
        ('CHEBI', 'CHEBI:3696')
    assert bioregistry.get_ns_id_from_bioregistry_curie('hgnc:1097') == \
        ('HGNC', '1097')


def test_get_bioregistry_prefix():
    assert bioregistry.get_bioregistry_prefix('PUBCHEM') == 'pubchem.compound'
    assert bioregistry.get_bioregistry_prefix('NXPFA') == 'nextprot.family'
    assert bioregistry.get_bioregistry_prefix('HGNC') == 'hgnc'


def test_get_bioregistry_curie():
    assert bioregistry.get_bioregistry_curie('PUBCHEM', '100101') == \
        'pubchem.compound:100101'
    assert bioregistry.get_bioregistry_curie('NXPFA', '01405') == \
        'nextprot.family:01405'
    assert bioregistry.get_bioregistry_curie('HGNC', '1097') == 'hgnc:1097'


def test_get_bioregistry_url():
    assert bioregistry.get_bioregistry_url('PUBCHEM', '100101') == \
        'https://bioregistry.io/pubchem.compound:100101'