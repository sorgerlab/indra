import os
from indra.preassembler.hierarchy_manager import HierarchyManager

def test_isa():
    rdf_file = os.path.join(os.path.dirname(__file__), 
        'entity_hierarchy.rdf')
    print rdf_file
    hm = HierarchyManager(rdf_file)
    assert(hm.isa('BRAF', 'RAF'))
