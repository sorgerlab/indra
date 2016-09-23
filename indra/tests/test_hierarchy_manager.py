from __future__ import print_function, unicode_literals
import os
from indra.preassembler.hierarchy_manager import hierarchies
from indra.statements import get_valid_location, InvalidLocationError

ent_hierarchy = hierarchies['entity']
mod_hierarchy = hierarchies['modification']
act_hierarchy = hierarchies['activity']
comp_hierarchy = hierarchies['cellular_component']

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
    assert comp_hierarchy.partof('INDRA', 'cytoplasm', 'INDRA', 'cell')

def test_partof_comp_not():
    assert not comp_hierarchy.partof('INDRA', 'cell', 'INDRA', 'cytoplasm')

def test_partof_comp_none():
    assert comp_hierarchy.partof('INDRA', 'cytoplasm', 'INDRA', None)

def test_partof_comp_none_none():
    assert comp_hierarchy.partof('INDRA', None, 'INDRA', None)

def test_partof_comp_none_not():
    assert not comp_hierarchy.partof('INDRA', None, 'INDRA', 'cytoplasm')

