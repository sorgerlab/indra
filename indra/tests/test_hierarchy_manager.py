import os
from copy import deepcopy
from indra.util import unicode_strs
from indra.preassembler.hierarchy_manager import hierarchies, \
    HierarchyManager, get_bio_hierarchies, YamlHierarchyManager, \
    get_wm_hierarchies
from indra.preassembler.make_wm_ontologies import eidos_ont_url, \
    rdf_graph_from_yaml, load_yaml_from_url


ent_hierarchy = hierarchies['entity']
mod_hierarchy = hierarchies['modification']
act_hierarchy = hierarchies['activity']
comp_hierarchy = hierarchies['cellular_component']
eidos_ns = 'https://github.com/clulab/eidos/wiki/JSON-LD/Grounding#'


def test_hierarchy_unicode():
    # Test all the hierarchies except the comp_hierarchy, which is an
    # RDF graph
    assert unicode_strs((ent_hierarchy.isa_closure,
                         ent_hierarchy.partof_closure))
    assert unicode_strs((mod_hierarchy.isa_closure,
                         mod_hierarchy.partof_closure))
    assert unicode_strs((act_hierarchy.isa_closure,
                         act_hierarchy.partof_closure))


def test_isa_entity():
    assert ent_hierarchy.isa('HGNC', '1097', 'FPLX', 'RAF')


def test_isa_entity2():
    assert not ent_hierarchy.isa('HGNC', '1097', 'HGNC', 'ARAF')


def test_isa_entity3():
    assert not ent_hierarchy.isa('FPLX', 'RAF', 'HGNC', '1097')


def test_partof_entity():
    assert ent_hierarchy.partof('FPLX', 'HIF_alpha', 'FPLX', 'HIF')


def test_isa_or_partof_entity():
    assert ent_hierarchy.isa_or_partof('HGNC', '9385', 'FPLX', 'AMPK')


def test_partof_entity_not():
    assert not ent_hierarchy.partof('FPLX', 'HIF1', 'FPLX', 'HIF_alpha')


def test_isa_mod():
    assert mod_hierarchy.isa('INDRA_MODS', 'phosphorylation',
                             'INDRA_MODS', 'modification')


def test_isa_mod_not():
    assert not mod_hierarchy.isa('INDRA_MODS', 'phosphorylation',
                                 'INDRA_MODS', 'ubiquitination')


def test_isa_activity():
    assert act_hierarchy.isa('INDRA_ACTIVITIES', 'kinase',
                             'INDRA_ACTIVITIES', 'activity')


def test_isa_activity_not():
    assert not act_hierarchy.isa('INDRA_ACTIVITIES', 'kinase',
                                 'INDRA_ACTIVITIES', 'phosphatase')


def test_partof_comp():
    assert comp_hierarchy.partof('INDRA_LOCATIONS', 'cytoplasm',
                                 'INDRA_LOCATIONS', 'cell')


def test_partof_comp_not():
    assert not comp_hierarchy.partof('INDRA_LOCATIONS', 'cell',
                                     'INDRA_LOCATIONS', 'cytoplasm')


def test_partof_comp_none():
    assert comp_hierarchy.partof('INDRA_LOCATIONS', 'cytoplasm',
                                 'INDRA_LOCATIONS', None)


def test_partof_comp_none_none():
    assert comp_hierarchy.partof('INDRA_LOCATIONS', None,
                                 'INDRA_LOCATIONS', None)


def test_partof_comp_none_not():
    assert not comp_hierarchy.partof('INDRA_LOCATIONS', None,
                                     'INDRA_LOCATIONS', 'cytoplasm')


def test_get_children():
    raf = 'http://identifiers.org/fplx/RAF'
    braf = 'http://identifiers.org/hgnc.symbol/BRAF'
    mapk = 'http://identifiers.org/fplx/MAPK'
    ampk = 'http://identifiers.org/fplx/AMPK'
    # Look up RAF
    rafs = ent_hierarchy.get_children(raf)
    # Should get three family members
    assert isinstance(rafs, list), rafs
    assert len(rafs) == 3
    assert unicode_strs(rafs)
    # The lookup of a gene-level entity should not return any additional
    # entities
    brafs = ent_hierarchy.get_children(braf)
    assert isinstance(brafs, list)
    assert len(brafs) == 0
    assert unicode_strs(brafs)
    mapks = ent_hierarchy.get_children(mapk)
    assert len(mapks) == 12, mapks
    assert unicode_strs(mapks)
    # Make sure we can also do this in a case involving both family and complex
    # relationships
    ampks = ent_hierarchy.get_children(ampk)
    assert len(ampks) == 22, ampks
    ag_none = ''
    none_children = ent_hierarchy.get_children('')
    assert isinstance(none_children, list)
    assert len(none_children) == 0


def test_mtorc_children():
    mtorc1 = 'http://identifiers.org/fplx/mTORC1'
    mtorc2 = 'http://identifiers.org/fplx/mTORC2'
    ch1 = ent_hierarchy.get_children(mtorc1)
    ch2 = ent_hierarchy.get_children(mtorc2)
    assert 'http://identifiers.org/hgnc.symbol/RICTOR' not in ch1
    assert 'http://identifiers.org/hgnc.symbol/RPTOR' not in ch2


def test_mtorc_get_parents():
    rictor = 'http://identifiers.org/hgnc/28611'  # RICTOR
    p = ent_hierarchy.get_parents(rictor, 'all')
    assert len(p) == 1
    assert list(p)[0] == 'http://identifiers.org/fplx/mTORC2'


def test_mtorc_transitive_closure():
    rictor = 'http://identifiers.org/hgnc/28611'  # RICTOR
    mtorc2 = 'http://identifiers.org/fplx/mTORC2'
    assert (rictor, mtorc2) in ent_hierarchy.partof_closure


def test_mtorc_partof_no_tc():
    ent_hierarchy_no_tc = deepcopy(ent_hierarchy)
    ent_hierarchy_no_tc.isa_closure = {}
    ent_hierarchy_no_tc.partof_closure = {}
    assert ent_hierarchy_no_tc.partof('HGNC', '30287', 'FPLX', 'mTORC1')
    assert not ent_hierarchy_no_tc.partof('HGNC', '30287', 'FPLX', 'mTORC2')


def test_erk_isa_no_tc():
    ent_hierarchy_no_tc = deepcopy(ent_hierarchy)
    ent_hierarchy_no_tc.isa_closure = {}
    ent_hierarchy_no_tc.partof_closure = {}
    assert ent_hierarchy_no_tc.isa('HGNC', '6871', 'FPLX', 'MAPK')
    assert not ent_hierarchy_no_tc.isa('HGNC', '6871', 'FPLX', 'JNK')


def test_get_parents():
    prkaa1 = 'http://identifiers.org/hgnc/9376'  # PRKAA1
    ampk = 'http://identifiers.org/fplx/AMPK'
    p1 = ent_hierarchy.get_parents(prkaa1, 'all')
    assert len(p1) == 8, p1
    assert ampk in p1
    p2 = ent_hierarchy.get_parents(prkaa1, 'immediate')
    assert len(p2) == 7, p2
    # This is to make sure we're getting an URI string
    assert unicode_strs(p2)
    assert ampk not in p2
    p3 = ent_hierarchy.get_parents(prkaa1, 'top')
    assert len(p3) == 1, p3
    assert ampk in p3


def test_chebi_isa():
    assert ent_hierarchy.isa('CHEBI', 'CHEBI:87307', 'CHEBI', 'CHEBI:36962')


def test_load_eid_hierarchy():
    eidos_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/eidos/eidos_ontology.rdf')
    hm = HierarchyManager(eidos_ont, True, True)
    assert hm.isa_closure
    eidos_isa = lambda a, b: hm.isa('UN', a, 'UN', b)
    assert eidos_isa('UN/events/human/conflict',
                     'UN/events/human')
    assert not eidos_isa('UN/events/human/conflict',
                         'UN/events/human/human_migration')
    assert eidos_isa('UN/entities/human/infrastructure',
                     'UN/entities')
    assert eidos_isa('UN/events/natural_disaster/storm',
                     'UN/events')
    assert not eidos_isa('UN/events',
                         'UN/events/natural/weather/storm')
    # Test case where graph is not given
    hm = HierarchyManager(None, True, True)
    hm.load_from_rdf_file(eidos_ont)
    assert eidos_isa('UN/events/natural_disaster/storm',
                     'UN/events')
    # Test loading from string
    with open(eidos_ont, 'r') as fh:
        hm = HierarchyManager(None, True, True)
        hm.load_from_rdf_string(fh.read())
    assert eidos_isa('UN/events/natural_disaster/storm',
                     'UN/events')
    # Test loading from Graph
    import rdflib
    g = rdflib.Graph()
    g.parse(eidos_ont, format='nt')
    hm = HierarchyManager(None, True, True)
    hm.load_from_rdf_graph(g)
    assert eidos_isa('UN/events/natural_disaster/storm',
                     'UN/events')


def test_load_trips_hierarchy():
    trips_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/cwms/trips_ontology.rdf')
    hm = HierarchyManager(trips_ont, True, True)
    assert hm.isa_closure
    trips_isa = lambda a, b: hm.isa('CWMS', a, 'CWMS', b)
    assert trips_isa('ONT::TRUCK', 'ONT::VEHICLE')
    assert not trips_isa('ONT::VEHICLE', 'ONT::TRUCK')
    assert trips_isa('ONT::MONEY', 'ONT::PHYS-OBJECT')
    assert trips_isa('ONT::TABLE', 'ONT::MANUFACTURED-OBJECT')


def test_load_sofia_hierarchy():
    sofia_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/sofia/sofia_ontology.rdf')
    hm = HierarchyManager(sofia_ont, True, True)
    assert hm.isa_closure
    sofia_isa = lambda a, b: hm.isa('SOFIA', a, 'SOFIA', b)
    assert sofia_isa('Accessibility/Accessibility', 'Accessibility')
    assert not sofia_isa('Movement/Transportation', 'Movement/Human_Migration')
    assert sofia_isa('Movement/Human_Migration', 'Movement')


def test_load_hume_hierarchy():
    hume_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            '../sources/hume/hume_ontology.rdf')
    hm = HierarchyManager(hume_ont, True, True)
    assert hm.isa_closure
    hume_isa = lambda a, b: hm.isa('HUME', a, 'HUME', b)
    assert hume_isa('entity/academic_discipline', 'entity')
    assert not hume_isa('entity', 'entity/academic_discipline')
    assert hume_isa('event/healthcare/human_disease',
                    'event/healthcare')


def test_same_components():
    uri_prkag1 = ent_hierarchy.get_uri('HGNC', '9385')  # PRKAG1
    uri_ampk = ent_hierarchy.get_uri('FPLX', 'AMPK')

    c1 = ent_hierarchy.components[uri_prkag1]
    c2 = ent_hierarchy.components[uri_ampk]
    assert c1 == c2


def test_bio_hierarchy_pickles():
    h1 = get_bio_hierarchies()
    h2 = get_bio_hierarchies(from_pickle=False)
    for key in h1.keys():
        assert len(h1[key].graph) == len(h2[key].graph)


def test_yaml_hm():
    yml = load_yaml_from_url(eidos_ont_url)
    hm = YamlHierarchyManager(yml, rdf_graph_from_yaml, True)

    entry = 'UN/events/natural_disaster/snowpocalypse'
    hm.add_entry(entry)
    assert hm.isa('UN', entry, 'UN', '/'.join(entry.split('/')[:-1]))

    entry = 'UN/events/galactic/alien_invasion'
    hm.add_entry(entry)
    assert hm.isa('UN', entry, 'UN', '/'.join(entry.split('/')[:-1]))
    assert hm.isa('UN', entry, 'UN', '/'.join(entry.split('/')[:-2]))


def test_hm_opposite():
    hierarchies = get_wm_hierarchies()
    concept1 = 'wm/concept/causal_factor/access/food_shortage'
    concept2 = ('wm/concept/causal_factor/economic_and_commerce/'
                'economic_activity/market/supply/food_supply')
    concept3 = ('wm/concept/causal_factor/environmental/meteorologic/'
                'precipitation/flooding')
    assert hierarchies['entity'].is_opposite('WM', concept1, 'WM', concept2)
    assert hierarchies['entity'].is_opposite('WM', concept2, 'WM', concept1)
    assert not hierarchies['entity'].is_opposite('WM', concept1, 'WM',
                                                 concept3)


def test_hm_equal():
    hierarchies = get_wm_hierarchies()
    concept1 = ('wm/concept/causal_factor/crisis_and_disaster/environmental/'
                'natural_disaster/flooding')
    concept2 = ('wm/concept/causal_factor/environmental/meteorologic/'
                'precipitation/flooding')
    concept3 = 'wm/concept/causal_factor/access/food_shortage'

    assert hierarchies['entity'].is_equal('WM', concept1, 'WM', concept2)
    assert hierarchies['entity'].is_equal('WM', concept2, 'WM', concept1)
    assert not hierarchies['entity'].is_equal('WM', concept1, 'WM', concept3)
