from indra.statements import Agent
from indra.ontology.bio import bio_ontology
from indra.databases import go_client, hgnc_client
from indra.ontology.standardize import \
    standardize_agent_name, standardize_db_refs, standardize_name_db_refs


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


def test_name_standardize_hgnc_up():
    a1 = Agent('x', db_refs={'HGNC': '9387'})
    standardize_agent_name(a1, True)
    assert a1.name == 'PRKAG3'
    a1 = Agent('x', db_refs={'UP': 'Q9UGI9'})
    standardize_agent_name(a1, True)
    assert a1.name == 'PRKAG3'
    a1 = Agent('x', db_refs={'UP': 'Q8BGM7'})
    standardize_agent_name(a1, True)
    assert a1.name == 'Prkag3'


def test_name_standardize_chebi():
    a1 = Agent('x', db_refs={'CHEBI': 'CHEBI:15996'})
    standardize_agent_name(a1, False)
    assert a1.name == 'GTP'


def test_name_standardize_go():
    a1 = Agent('x', db_refs={'GO': 'GO:0006915'})
    standardize_agent_name(a1, False)
    assert a1.name == 'apoptotic process'


def test_name_standardize_mesh():
    a1 = Agent('x', db_refs={'MESH': 'D008545'})
    standardize_agent_name(a1, False)
    assert a1.name == 'Melanoma', a1.name


def test_name_standardize_mesh_go():
    a1 = Agent('x', db_refs={'MESH': 'D058750'})
    standardize_agent_name(a1, True)
    assert a1.db_refs['GO'] == 'GO:0001837'
    assert a1.name == 'epithelial to mesenchymal transition', a1.name
    a1 = Agent('x', db_refs={'GO': 'GO:0001837'})
    standardize_agent_name(a1, True)
    assert a1.db_refs['MESH'] == 'D058750'
    assert a1.name == 'epithelial to mesenchymal transition', a1.name


def test_name_standardize_mesh_other_db():
    a1 = Agent('x', db_refs={'MESH': 'D001194'})
    standardize_agent_name(a1, True)
    assert a1.db_refs['CHEBI'] == 'CHEBI:46661'
    assert a1.name == 'asbestos', a1.name

    db_refs = {'MESH': 'D000067777'}
    db_refs = standardize_db_refs(db_refs)
    assert db_refs.get('HGNC') == '3313', db_refs
    assert db_refs.get('UP') == 'Q12926', db_refs
    a2 = Agent('x', db_refs=db_refs)
    standardize_agent_name(a2)
    assert a2.name == 'ELAVL2'


def test_standardize_db_refs_efo_hp_doid():
    refs = standardize_db_refs({'EFO': '0009502'})
    assert refs.get('MESH') == 'D000007', refs
    refs = standardize_db_refs({'MESH': 'D000007'})
    assert refs.get('EFO') == '0009502', refs

    refs = standardize_db_refs({'HP': 'HP:0031801'})
    assert refs.get('MESH') == 'D064706', refs
    refs = standardize_db_refs({'MESH': 'D064706'})
    assert refs.get('HP') == 'HP:0031801', refs

    # Currently there is no one-to-many mapping in the direction towards MeSH
    # (there used to be) if there is again, we should test it here
    #refs = standardize_db_refs({'DOID': 'DOID:0060695'})
    #assert 'MESH' not in refs

    # One-to-many mappings away from MESH
    refs = standardize_db_refs({'MESH': 'D000071017'})
    assert 'DOID' not in refs

    refs = standardize_db_refs({'DOID': 'DOID:0060495'})
    assert refs.get('MESH') == 'D000067208'

    # This is an xrefs-based mapping that isn't in Gilda's resource file
    refs = standardize_db_refs({'EFO': '0000694'})
    assert refs.get('MESH') == 'D045169'


def test_standardize_name_efo_hp_doid():
    ag = Agent('x', db_refs={'HP': 'HP:0031801'})
    standardize_agent_name(ag)
    # Name based on MESH mapping
    assert ag.name == 'Vocal Cord Dysfunction'

    ag = Agent('x', db_refs={'HP': 'HP:0000002'})
    standardize_agent_name(ag)
    # Name based on HP itself
    assert ag.name == 'Abnormality of body height'

    ag = Agent('x', db_refs={'DOID': 'DOID:0014667'})
    standardize_agent_name(ag)
    # Name based on MESH mapping
    assert ag.name == 'Metabolic Diseases'

    ag = Agent('x', db_refs={'EFO': '1002050'})
    standardize_agent_name(ag)
    # Name based on MESH mapping
    assert ag.name == 'Nephritis', (ag.name, ag.db_refs)

    ag = Agent('x', db_refs={'EFO': '0000001'})
    standardize_agent_name(ag)
    # Name based on EFO itself
    assert ag.name == 'experimental factor', (ag.name, ag.db_refs)


def test_standardize_uppro():
    ag = Agent('x', db_refs={'UP': 'P01019'})
    standardize_agent_name(ag)
    assert ag.name == 'AGT'
    ag = Agent('x', db_refs={'UPPRO': 'PRO_0000032458'})
    standardize_agent_name(ag)
    assert ag.name == 'Angiotensin-2', ag.name
    ag = Agent('x', db_refs={'UPPRO': 'PRO_0000032458', 'UP': 'P01019'})
    standardize_agent_name(ag)
    assert ag.name == 'Angiotensin-2', ag.name


def test_uppro_fallback():
    # This UP chain has no name currently so we can test that the fallback
    # to naming by the UP ID is working
    ag = Agent('x', db_refs={'UP': 'Q6IE75', 'UPPRO': 'PRO_0000383648'})
    standardize_agent_name(ag)
    assert ag.name == 'Bace2'


def test_mirna_standardize():
    name, db_refs = standardize_name_db_refs({'HGNC': '31476'})
    assert db_refs['HGNC'] == '31476'
    assert db_refs['MIRBASE'] == 'MI0000060'
    assert name == 'MIRLET7A1'

    name, db_refs = standardize_name_db_refs({'MIRBASE': 'MI0001730'})
    assert db_refs['MIRBASE'] == 'MI0001730'
    assert name == 'mmu-mir-451a'


def test_drugbank_mappings():
    name, db_refs = standardize_name_db_refs({'DRUGBANK': 'DB00001'})
    assert db_refs.get('CHEBI') == 'CHEBI:142437', db_refs
    assert db_refs.get('CHEMBL') == 'CHEMBL1201666', db_refs
    assert name == 'lepirudin'
    # Here we test for alternative prioritization of name spaces
    name, db_refs = standardize_name_db_refs({'DRUGBANK': 'DB00001'},
                                             ns_order=['DRUGBANK', 'CHEBI'])
    # We expect to get the Drugbank standard name
    assert name == 'Lepirudin'


def test_standardize_up_isoform():
    refs = standardize_db_refs({'UP': 'Q99490'})
    assert refs == {'UP': 'Q99490', 'HGNC': '16921',
                    'EGID': '116986', 'MESH': 'C485997'}, refs
    refs = standardize_db_refs({'UP': 'Q99490-123'})
    assert refs == {'UP': 'Q99490-123', 'HGNC': '16921',
                    'EGID': '116986', 'MESH': 'C485997'}, refs


def test_standardize_chembl():
    db_refs = standardize_db_refs({'DRUGBANK': 'DB00305'})
    assert 'CHEMBL' in db_refs, db_refs
    assert db_refs['CHEMBL'] == 'CHEMBL105', db_refs


def test_efo_bfo_relations():
    assert set(bio_ontology.get_parents('EFO', '0004542')) == \
        {('BFO', '0000015'), ('EFO', '0000001')}


def test_name_lookup_obsolete():
    # This is a regression test to make sure we don't return another node
    # with the same name but which is obsolete (HGNC:11093)
    assert bio_ontology.get_id_from_name('HGNC', 'ALDH3A2') == \
        ('HGNC', '403')


def test_chebi_refinements():
    assert bio_ontology.partof('CHEBI', 'CHEBI:136692',
                               'CHEBI', 'CHEBI:365')
    assert not bio_ontology.partof('CHEBI', 'CHEBI:365',
                                   'CHEBI', 'CHEBI:136692')


def test_standardize_hgnc_fplx_mesh_bug():
    refs = standardize_db_refs({'HGNC': '1514'})
    assert refs['UP'] == 'P41180'
    assert refs['MESH'] == 'C095550'
    assert 'FPLX' not in refs

    refs = standardize_db_refs({'FPLX': 'Calcium_sensing_receptors'})
    assert refs['MESH'] == 'D044169'
    assert refs['HGNC_GROUP'] == '279'
    assert 'HGNC' not in refs


def test_ido_parents():
    parents = bio_ontology.get_parents('IDO', '0000514')
    assert ('IDO', '0000509') in parents


def test_lspci():
    assert bio_ontology.get_name('LSPCI', '18') == 'Pentane-1,5-Diamine'
    members = bio_ontology.get_children('LSPCI', '18')
    # These are some of the members, not all
    expected_members = {('CAS', '462-94-2'),
                        ('CHEBI', 'CHEBI:18127'),
                        ('CHEMBL', 'CHEMBL119296'),
                        ('PUBCHEM', '273')}
    assert expected_members < set(members)
