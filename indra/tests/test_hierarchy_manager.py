import os
from indra.preassembler.hierarchy_manager import HierarchyManager

rdf_file = os.path.join(os.path.dirname(__file__), 
        '../preassembler/entity_hierarchy.rdf')

def test_find_entity():
    hm = HierarchyManager(rdf_file)
    assert(hm.find_entity('BRAF'))

def test_find_entity2():
    hm = HierarchyManager(rdf_file)
    assert(hm.find_entity('abcdefghxyz') is None)

def test_find_entity3():
    hm = HierarchyManager(rdf_file)
    assert(hm.find_entity('RAF_FAMILY'))

def test_find_entity4():
    hm = HierarchyManager(rdf_file)
    assert(hm.find_entity('B-RAF1'))

def test_isa():
    hm = HierarchyManager(rdf_file)
    assert(hm.isa('BRAF', 'RAF'))

def test_isa2():
    hm = HierarchyManager(rdf_file)
    assert(hm.isa('B-RAF1', 'RAF_FAMILY'))

def test_isa3():
    hm = HierarchyManager(rdf_file)
    assert(not hm.isa('BRAF', 'ARAF'))

def test_isa4():
    hm = HierarchyManager(rdf_file)
    assert(not hm.isa('RAF', 'BRAF'))
