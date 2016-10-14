from indra.preassembler import make_activity_hierarchy as make_act
from indra.preassembler import make_entity_hierarchy as make_ent
from indra.preassembler import make_modification_hierarchy as make_mod
from indra.preassembler import make_cellular_component_hierarchy as make_comp

def test_make_activity_hierarchy():
    make_act.main()

def test_make_entity_hierarchy():
    make_ent.main(make_ent.relations_file)

def test_make_modification_hierarchy():
    make_mod.main()

"""
def test_make_cellular_component_hierarchy():
    make_comp.main()
"""
