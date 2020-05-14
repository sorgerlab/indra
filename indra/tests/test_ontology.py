from indra.ontology.bio import bio_ontology
from indra.ontology.world import world_ontology
from indra.databases import go_client, hgnc_client


def test_isa_entity():
    assert bio_ontology.isa('HGNC', '1097', 'FPLX', 'RAF')


def test_isa_entity2():
    assert not bio_ontology.isa('HGNC', '1097', 'HGNC', '646')


def test_isa_entity3():
    assert not bio_ontology.isa('FPLX', 'RAF', 'HGNC', '1097')


def test_partof_entity():
    assert bio_ontology.partof('FPLX', 'HIF_alpha', 'FPLX', 'HIF')


def test_isa_or_partof_entity():
    assert bio_ontology.isa_or_partof('HGNC', '9385', 'FPLX', 'AMPK')


def test_partof_entity_not():
    assert not bio_ontology.partof('FPLX', 'HIF1', 'FPLX', 'HIF_alpha')


def test_isa_mod():
    assert bio_ontology.isa('INDRA_MODS', 'phosphorylation',
                            'INDRA_MODS', 'modification')


def test_isa_mod_not():
    assert not bio_ontology.isa('INDRA_MODS', 'phosphorylation',
                                'INDRA_MODS', 'ubiquitination')


def test_isa_activity():
    assert bio_ontology.isa('INDRA_ACTIVITIES', 'kinase',
                            'INDRA_ACTIVITIES', 'activity')


def test_isa_activity_not():
    assert not bio_ontology.isa('INDRA_ACTIVITIES', 'kinase',
                                'INDRA_ACTIVITIES', 'phosphatase')


def test_partof_comp():
    assert bio_ontology.isa_or_partof(
        'GO', go_client.get_go_id_from_label('cytoplasm'),
        'GO', go_client.get_go_id_from_label('cellular_component'))


def test_partof_comp_not():
    assert not bio_ontology.isa_or_partof(
        'GO', go_client.get_go_id_from_label('cellular_component'),
        'GO', go_client.get_go_id_from_label('cytoplasm'))


def test_get_children():
    rafs = bio_ontology.get_children('FPLX', 'RAF')
    assert isinstance(rafs, list), rafs
    assert len(rafs) == 3, rafs

    brafs = bio_ontology.get_children('HGNC', hgnc_client.get_hgnc_id('BRAF'))
    assert isinstance(brafs, list), brafs
    assert len(brafs) == 0

    mapks = bio_ontology.get_children('FPLX', 'MAPK')
    assert len(mapks) == 12, mapks

    ampks = bio_ontology.get_children('FPLX', 'AMPK')
    assert len(ampks) == 22, ampks


def test_mtorc_children():
    ch1 = bio_ontology.get_children('FPLX', 'mTORC1')
    ch2 = bio_ontology.get_children('FPLX', 'mTORC2')
    assert ('HGNC', hgnc_client.get_hgnc_id('RICTOR')) not in ch1
    assert ('HGNC', hgnc_client.get_hgnc_id('RPTOR')) not in ch2


def test_mtorc_get_parents():
    p = bio_ontology.get_parents('HGNC', hgnc_client.get_hgnc_id('RICTOR'))
    assert len(p) == 1
    assert p == [('FPLX', 'mTORC2')]


def test_mtorc_transitive_closure():
    assert bio_ontology.partof('HGNC', hgnc_client.get_hgnc_id('RICTOR'),
                               'FPLX', 'mTORC2')
    assert not bio_ontology.partof('HGNC', hgnc_client.get_hgnc_id('RPTOR'),
                                   'FPLX', 'mTORC2')


def test_erk_isa():
    assert bio_ontology.isa('HGNC', '6871', 'FPLX', 'MAPK')
    assert not bio_ontology.isa('HGNC', '6871', 'FPLX', 'JNK')


def test_get_parents():
    prkaa1 = ('HGNC', '9376')
    ampk = ('FPLX', 'AMPK')
    p1 = bio_ontology.get_parents(*prkaa1)
    assert len(p1) == 8, p1
    assert ampk in p1

    # FIXME: implement these
    # p2 = ent_hierarchy.get_parents(prkaa1, 'immediate')
    # assert len(p2) == 7, p2
    # This is to make sure we're getting an URI string
    # assert unicode_strs(p2)
    # assert ampk not in p2
    # p3 = ent_hierarchy.get_parents(prkaa1, 'top')
    # assert len(p3) == 1, p3
    # assert ampk in p3


def test_chebi_isa():
    assert bio_ontology.isa('CHEBI', 'CHEBI:87307', 'CHEBI', 'CHEBI:36962')

# FIXME: implement components
# def test_same_components():
#    uri_prkag1 = ent_hierarchy.get_uri('HGNC', '9385')  # PRKAG1
#    uri_ampk = ent_hierarchy.get_uri('FPLX', 'AMPK')
#
#    c1 = ent_hierarchy.components[uri_prkag1]
#    c2 = ent_hierarchy.components[uri_ampk]
#    assert c1 == c2


def test_hm_opposite_polarity():
    concept1 = 'wm/concept/causal_factor/food_insecurity/food_instability'
    concept2 = 'wm/concept/causal_factor/food_security/food_stability'
    concept3 = ('wm/concept/causal_factor/environmental/meteorologic/'
                'precipitation/flooding')
    assert world_ontology.is_opposite('WM', concept1, 'WM', concept2)
    assert world_ontology.is_opposite('WM', concept2, 'WM', concept1)
    assert not world_ontology.is_opposite('WM', concept1, 'WM', concept3)
    assert world_ontology.get_polarity('WM', concept1) == -1
    assert world_ontology.get_polarity('WM', concept2) == 1
    assert world_ontology.get_polarity('UN', 'something') is None
