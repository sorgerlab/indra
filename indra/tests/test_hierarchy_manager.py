from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from copy import deepcopy
from indra.preassembler.hierarchy_manager import hierarchies
from indra.statements import get_valid_location, InvalidLocationError, Agent
from indra.util import unicode_strs

ent_hierarchy = hierarchies['entity']
mod_hierarchy = hierarchies['modification']
act_hierarchy = hierarchies['activity']
comp_hierarchy = hierarchies['cellular_component']

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
    assert(ent_hierarchy.isa('HGNC', 'BRAF', 'FPLX', 'RAF'))

def test_isa_entity2():
    assert(not ent_hierarchy.isa('HGNC', 'BRAF', 'HGNC', 'ARAF'))

def test_isa_entity3():
    assert(not ent_hierarchy.isa('FPLX', 'RAF', 'HGNC', 'BRAF'))

def test_partof_entity():
    assert ent_hierarchy.partof('FPLX', 'HIF_alpha', 'FPLX', 'HIF')

def test_partof_entity_not():
    assert not ent_hierarchy.partof('FPLX', 'HIF1', 'FPLX', 'HIF_alpha')

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

def test_get_children():
    raf = 'http://sorger.med.harvard.edu/indra/entities/RAF'
    braf = 'http://identifiers.org/hgnc.symbol/BRAF'
    mapk = 'http://sorger.med.harvard.edu/indra/entities/MAPK'
    ampk = 'http://sorger.med.harvard.edu/indra/entities/AMPK'
    # Look up RAF
    rafs = ent_hierarchy.get_children(raf)
    # Should get three family members
    assert isinstance(rafs, list)
    assert len(rafs) == 3
    assert unicode_strs(rafs)
    # The lookup of a gene-level entity should not return any additional
    # entities
    brafs = ent_hierarchy.get_children(braf)
    assert isinstance(brafs, list)
    assert len(brafs) == 0
    assert unicode_strs(brafs)
    mapks = ent_hierarchy.get_children(mapk)
    assert len(mapks) == 12
    assert unicode_strs(mapks)
    # Make sure we can also do this in a case involving both family and complex
    # relationships
    ampks = ent_hierarchy.get_children(ampk)
    assert len(ampks) == 22
    ag_none = ''
    none_children = ent_hierarchy.get_children('')
    assert isinstance(none_children, list)
    assert len(none_children) == 0

def test_mtorc_children():
    mtorc1 = 'http://sorger.med.harvard.edu/indra/entities/mTORC1'
    mtorc2 = 'http://sorger.med.harvard.edu/indra/entities/mTORC2'
    ch1 = ent_hierarchy.get_children(mtorc1)
    ch2 = ent_hierarchy.get_children(mtorc2)
    assert('http://identifiers.org/hgnc.symbol/RICTOR' not in ch1)
    assert('http://identifiers.org/hgnc.symbol/RPTOR' not in ch2)

def test_mtorc_get_parents():
    rictor = 'http://identifiers.org/hgnc.symbol/RICTOR'
    p = ent_hierarchy.get_parents(rictor, 'all')
    assert(len(p) == 1)
    assert(list(p)[0] == 'http://sorger.med.harvard.edu/indra/entities/mTORC2')

def test_mtorc_transitive_closure():
    rictor = 'http://identifiers.org/hgnc.symbol/RICTOR'
    p = ent_hierarchy.partof_closure.get(rictor)
    assert(len(p) == 1)
    assert(p[0] == 'http://sorger.med.harvard.edu/indra/entities/mTORC2')

def test_mtorc_partof_no_tc():
    ent_hierarchy_no_tc = deepcopy(ent_hierarchy)
    ent_hierarchy_no_tc.isa_closure = {}
    ent_hierarchy_no_tc.partof_closure = {}
    assert(ent_hierarchy_no_tc.partof('HGNC', 'RPTOR', 'FPLX', 'mTORC1'))
    assert(not ent_hierarchy_no_tc.partof('HGNC', 'RPTOR', 'FPLX', 'mTORC2'))

def test_erk_isa_no_tc():
    ent_hierarchy_no_tc = deepcopy(ent_hierarchy)
    ent_hierarchy_no_tc.isa_closure = {}
    ent_hierarchy_no_tc.partof_closure = {}
    assert(ent_hierarchy_no_tc.isa('HGNC', 'MAPK1', 'FPLX', 'MAPK'))
    assert(not ent_hierarchy_no_tc.isa('HGNC', 'MAPK1', 'FPLX', 'JNK'))

def test_get_parents():
    prkaa1 = 'http://identifiers.org/hgnc.symbol/PRKAA1'
    ampk = 'http://sorger.med.harvard.edu/indra/entities/AMPK'
    p1 = ent_hierarchy.get_parents(prkaa1, 'all')
    assert(len(p1) == 8)
    assert(ampk in p1)
    p2 = ent_hierarchy.get_parents(prkaa1, 'immediate')
    assert(len(p2) == 7)
    assert (ampk not in p2)
    p3 = ent_hierarchy.get_parents(prkaa1, 'top')
    assert(len(p3) == 1)
    assert (ampk in p3)
