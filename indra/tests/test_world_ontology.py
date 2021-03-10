import copy
from indra.ontology.world import world_ontology
from indra.ontology.world.ontology import WorldOntology


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


def test_world_ontology_add_entry():
    ont = copy.deepcopy(world_ontology)
    nat_dis = ('wm/concept/causal_factor/crisis_and_disaster/'
               'environmental_disasters/natural_disaster')

    new_node = nat_dis + '/floods'
    assert not ont.isa('WM', new_node, 'WM', nat_dis)
    ont.add_entry(new_node, examples=['floods'])
    assert ont.isa('WM', new_node, 'WM', nat_dis)
    ont_yml = ont.dump_yml_str()


def test_load_intermediate_nodes():
    url = ('https://raw.githubusercontent.com/clulab/eidos/posneg/src/main/'
           'resources/org/clulab/wm/eidos/english/ontologies/'
           'wm_posneg_metadata.yml')
    wo = WorldOntology(url)
    wo.initialize()
