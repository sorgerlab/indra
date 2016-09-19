import os
from indra.preassembler.hierarchy_manager import HierarchyManager
from indra.statements import get_valid_location, InvalidLocationError

entity_file = os.path.join(os.path.dirname(__file__), 
        '../resources/entity_hierarchy.rdf')

mod_file = os.path.join(os.path.dirname(__file__), 
        '../resources/modification_hierarchy.rdf')

act_file = os.path.join(os.path.dirname(__file__), 
        '../resources/activity_hierarchy.rdf')

comp_file = os.path.join(os.path.dirname(__file__), 
        '../resources/cellular_component_hierarchy.rdf')

ent_hierarchy = HierarchyManager(entity_file)
mod_hierarchy = HierarchyManager(mod_file)
act_hierarchy = HierarchyManager(act_file)
comp_hierarchy = HierarchyManager(comp_file)

def test_isa_entity():
    assert(ent_hierarchy.isa('HGNC', 'BRAF', 'BE', 'RAF'))

def test_isa_entity2():
    assert(not ent_hierarchy.isa('HGNC', 'BRAF', 'HGNC', 'ARAF'))

def test_isa_entity3():
    assert(not ent_hierarchy.isa('BE', 'RAF', 'HGNC', 'BRAF'))

def test_partof_entity():
    assert ent_hierarchy.partof('BE', 'HIF1_alpha', 'BE', 'HIF1')

def test_partof_entity_not():
    assert not ent_hierarchy.partof('BE', 'HIF1', 'BE', 'HIF1_alpha')

def test_isa_mod():
    assert(mod_hierarchy.isa('INDRA', 'phosphorylation',
                             'INDRA', 'modification'))

def test_isa_mod_not():
    assert(not mod_hierarchy.isa('INDRA', 'phosphorylation',
                                 'INDRA', 'ubiquitination'))

def test_isa_activity():
    assert act_hierarchy.isa('INDRA', 'kinase', 'INDRA', 'activity')

def test_isa_activity_not():
    assert not act_hierarchy.isa('INDRA', 'kinase', 'INDRA', 'phosphatase')

def test_partof_comp():
    cyto_loc = get_valid_location('cytoplasm')
    cell_loc = get_valid_location('cell')
    assert comp_hierarchy.partof('INDRA', cyto_loc, 'INDRA', cell_loc)

def test_partof_comp_go_id():
    cyto_loc = get_valid_location('GO:0005737')
    cell_loc = get_valid_location('GO:0005623')
    assert comp_hierarchy.partof('INDRA', cyto_loc, 'INDRA', cell_loc)

def test_partof_comp_not():
    cyto_loc = get_valid_location('cytoplasm')
    cell_loc = get_valid_location('cell')
    assert not comp_hierarchy.partof('INDRA', cell_loc, 'INDRA', cyto_loc)

def test_partof_comp_none():
    cyto_loc = get_valid_location('cytoplasm')
    assert comp_hierarchy.partof('INDRA', cyto_loc, 'INDRA', None)

def test_partof_comp_none_none():
    assert comp_hierarchy.partof('INDRA', None, 'INDRA', None)

def test_partof_comp_none_not():
    cyto_loc = get_valid_location('cytoplasm')
    assert not comp_hierarchy.partof('INDRA', None, 'INDRA', cyto_loc)

