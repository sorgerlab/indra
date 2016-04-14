import os
from indra.preassembler.hierarchy_manager import HierarchyManager

entity_file = os.path.join(os.path.dirname(__file__), 
        '../resources/entity_hierarchy.rdf')

mod_file = os.path.join(os.path.dirname(__file__), 
        '../resources/modification_hierarchy.rdf')

def test_find_entity():
    hm = HierarchyManager(entity_file)
    assert(hm.find_entity('BRAF'))

def test_find_entity2():
    hm = HierarchyManager(entity_file)
    assert(hm.find_entity('abcdefghxyz') is None)

def test_find_entity3():
    hm = HierarchyManager(entity_file)
    assert(hm.find_entity('RAF_FAMILY'))

def test_find_entity4():
    hm = HierarchyManager(entity_file)
    assert(hm.find_entity('B-RAF1'))

def test_find_mod():
    hm = HierarchyManager(mod_file)
    assert(hm.find_entity('phosphorylation'))

def test_find_mod2():
    hm = HierarchyManager(mod_file)
    assert(hm.find_entity('sumoylation'))

def test_isa():
    hm = HierarchyManager(entity_file)
    assert(hm.isa('BRAF', 'RAF'))

def test_isa2():
    hm = HierarchyManager(entity_file)
    assert(hm.isa('B-RAF1', 'RAF_FAMILY'))

def test_isa3():
    hm = HierarchyManager(entity_file)
    assert(not hm.isa('BRAF', 'ARAF'))

def test_isa4():
    hm = HierarchyManager(entity_file)
    assert(not hm.isa('RAF', 'BRAF'))

def test_isa_mod():
    hm = HierarchyManager(mod_file)
    assert(hm.isa('phosphorylation', 'modification'))

def test_isa_mod_not():
    hm = HierarchyManager(mod_file)
    assert(not hm.isa('phosphorylation', 'ubiquitination'))
