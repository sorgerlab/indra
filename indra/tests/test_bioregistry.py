from indra.databases import bioregistry_client


def test_get_ns_from_bioregistry():
    assert bioregistry_client.get_ns_from_bioregistry('xxxx') is None
    assert bioregistry_client.get_ns_from_bioregistry('noncodev4.rna') == \
        'NONCODE'
    assert bioregistry_client.get_ns_from_bioregistry('chebi') == 'CHEBI'


def test_get_ns_id_from_bioregistry():
    assert bioregistry_client.get_ns_id_from_bioregistry('xxxx', 'xxxx') == \
        (None, None)
    assert bioregistry_client.get_ns_id_from_bioregistry('chebi', '3696') == \
        ('CHEBI', 'CHEBI:3696')
    assert bioregistry_client.get_ns_id_from_bioregistry('hgnc', '1097') == \
        ('HGNC', '1097')
    assert bioregistry_client.get_ns_id_from_bioregistry('cellosaurus',
                                                         '1234') == \
        ('CVCL', 'CVCL_1234')


def test_get_ns_id_from_bioregistry_curie():
    assert bioregistry_client.get_ns_id_from_bioregistry_curie('xxxx:xxxx') == \
        (None, None)
    assert bioregistry_client.get_ns_id_from_bioregistry_curie('chebi:3696') == \
        ('CHEBI', 'CHEBI:3696')
    assert bioregistry_client.get_ns_id_from_bioregistry_curie('hgnc:1097') == \
        ('HGNC', '1097')


def test_get_bioregistry_prefix():
    assert bioregistry_client.get_bioregistry_prefix('PUBCHEM') == \
        'pubchem.compound'
    assert bioregistry_client.get_bioregistry_prefix('NXPFA') == \
        'nextprot.family'
    assert bioregistry_client.get_bioregistry_prefix('HGNC') == 'hgnc'


def test_get_bioregistry_curie():
    assert bioregistry_client.get_bioregistry_curie('PUBCHEM', '100101') == \
        'pubchem.compound:100101'
    assert bioregistry_client.get_bioregistry_curie('NXPFA', '01405') == \
        'nextprot.family:01405'
    assert bioregistry_client.get_bioregistry_curie('HGNC', '1097') == \
        'hgnc:1097'
    assert bioregistry_client.get_bioregistry_curie('CVCL', 'CVCL_1234') == \
        'cellosaurus:1234'


def test_get_bioregistry_url():
    assert bioregistry_client.get_bioregistry_url('PUBCHEM', '100101') == \
        'https://bioregistry.io/pubchem.compound:100101'


def test_ensure_prefix_if_needed():
    assert bioregistry_client.ensure_prefix_if_needed('PUBCHEM', '100101') == \
        '100101'
    assert bioregistry_client.ensure_prefix_if_needed('CHEBI', '3696') == \
        'CHEBI:3696'
    assert bioregistry_client.ensure_prefix_if_needed('CHEBI', 'CHEBI:3696') == \
        'CHEBI:3696'
