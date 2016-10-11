from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from nose.tools import raises
from indra.preassembler.hierarchy_manager import HierarchyManager
from indra.preassembler.hierarchy_manager import hierarchies
from indra.statements import *
from indra.util import unicode_strs

# Argument checking for ActiveForms ----------------------------

def test_activitymod_sitelist_of_ints():
    """Check that mod positions specified as ints get promoted to strings."""
    st = ActiveForm(Agent('MAP2K1', mods=
                              [ModCondition('phosphorylation', 'serine', 218),
                               ModCondition('phosphorylation', 'serine', 222)]),
                              'kinase', True)
    assert not isinstance(st.agent.mods[0].position, int)
    assert not isinstance(st.agent.mods[1].position, int)
    assert unicode_strs(st)

def testactivitymod_string_string():
    """Check that string mod position is preserved"""
    st = ActiveForm(Agent('MAP2K1', mods=
                            [ModCondition('phosphorylation', 'serine', '218')]),
                             'kinase', True)
    assert not isinstance(st.agent.mods[0].position, int)
    assert unicode_strs(st)

def testactivitymod_string_none():
    """Check that None mod position is preserved"""
    st = ActiveForm(Agent('MAP2K1', mods=
                          [ModCondition('phosphorylation', 'serine', None)]),
                          'kinase', True)
    assert (st.agent.mods[0].position is None)
    assert unicode_strs(st)

def testactivitymod_nolist():
    """Make sure mod is correctly turned into a list if it's
    a single ModCondition"""
    mc = ModCondition('phosphorylation')
    st = ActiveForm(Agent('MAP2K1', mods=mc),
                          'kinase', True)
    assert isinstance(st.agent.mods, list)
    assert unicode_strs(st)
    assert unicode_strs(mc)

# Checking for exact matching (except Evidence) between Agents/stmts ---------

def test_matches():
    ras = Agent('Ras')
    raf = Agent('Raf')
    st1 = Complex([ras, raf])
    st2 = Complex([ras, raf])
    assert(st1.matches(st2))
    assert unicode_strs(st1)

def test_matches_key():
    ras = Agent('Ras')
    raf = Agent('Raf')
    st1 = Complex([ras, raf])
    st2 = Complex([ras, raf])
    assert(st1.matches_key() == st2.matches_key())
    assert unicode_strs(st1)

def test_matches_key_unicode():
    ras = Agent('Ras')
    rasu = Agent(u'Ras')
    raf = Agent('Raf')
    st1 = Complex([ras, raf])
    st2 = Complex([rasu, raf])
    assert(st1.matches_key() == st2.matches_key())
    assert unicode_strs(st1)
    assert unicode_strs(st2)

def test_matches_key_unicode2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek, u'S')
    st2 = Phosphorylation(raf, mek, 'S')
    assert(st1.matches_key() == st2.matches_key())
    assert unicode_strs(st1)
    assert unicode_strs(st2)

def test_matches_key_unicode3():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek, 'S', u'222')
    st2 = Phosphorylation(raf, mek, 'S', '222')
    assert(st1.matches_key() == st2.matches_key())
    assert unicode_strs(st1)
    assert unicode_strs(st2)

def test_matches2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek)
    assert(st1.matches(st2))
    assert unicode_strs(st1)
    assert unicode_strs(st2)

def test_matches_key2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek)
    assert(st1.matches_key() == st2.matches_key())
    assert unicode_strs(st1)
    assert unicode_strs(st2)

def test_not_matches():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek, 'tyrosine')
    assert(not st1.matches(st2))
    assert unicode_strs(st1)
    assert unicode_strs(st2)

def test_not_matches_key():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek, 'tyrosine')
    assert(st1.matches_key() != st2.matches_key())
    assert unicode_strs(st1)
    assert unicode_strs(st2)

def test_matches_dbrefs():
    hras1 = Agent('HRAS', db_refs={'hgnc': 1111})
    hras2 = Agent('HRAS', db_refs={'hgnc': 9999})
    assert(hras1.matches(hras2))
    assert unicode_strs(hras1)
    assert unicode_strs(hras2)

def test_matches_key_dbrefs():
    hras1 = Agent('HRAS', db_refs={'hgnc': 1111})
    hras2 = Agent('HRAS', db_refs={'hgnc': 9999})
    assert(hras1.matches_key() == hras2.matches_key())
    assert unicode_strs((hras1, hras2))

def test_matches_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    assert(hras1.matches(hras2))
    assert unicode_strs((hras1, hras2))

def test_matches_key_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    assert(hras1.matches_key() == hras2.matches_key())
    assert unicode_strs((hras1, hras2))

def test_not_matches_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('RAF1'), True)])
    assert(not hras1.matches(hras2))
    assert unicode_strs((hras1, hras2))

def test_not_matches_key_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches_key() != hras2.matches_key())
    assert unicode_strs((hras1, hras2))

def test_not_matches_bound2():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), False)])
    assert(not hras1.matches(hras2))
    assert unicode_strs((hras1, hras2))

def test_not_matches_key_bound2():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), False)])
    assert(hras1.matches_key() != hras2.matches_key())
    assert unicode_strs((hras1, hras2))

def test_matches_bound_multiple():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                           BoundCondition(Agent('RAF1'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                           BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches(hras2))
    assert unicode_strs((hras1, hras2))

def test_matches_key_bound_multiple():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                           BoundCondition(Agent('RAF1'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                           BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches_key() == hras2.matches_key())
    assert unicode_strs((hras1, hras2))

def test_matches_bound_multiple_order():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('RAF1'), True),
                                           BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                           BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches(hras2))
    assert unicode_strs((hras1, hras2))

def test_matches_key_bound_multiple_order():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('RAF1'), True),
                                           BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                           BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches_key() == hras2.matches_key())
    assert unicode_strs((hras1, hras2))

def test_matches_agent_mod_order():
    hras1 = Agent('MAP2K1',
        mods=[ModCondition('phosphorylation'), ModCondition('ubiquitination')])
    hras2 = Agent('MAP2K1',
        mods=[ModCondition('ubiquitination'), ModCondition('phosphorylation')])
    assert(hras1.matches(hras2))
    assert unicode_strs((hras1, hras2))

def test_refinement_agent_mod_order():
    hras1 = Agent('MAP2K1',
        mods=[ModCondition('phosphorylation', 'S'),
              ModCondition('ubiquitination')])
    hras2 = Agent('MAP2K1',
        mods=[ModCondition('ubiquitination'), ModCondition('phosphorylation')])
    assert(hras1.refinement_of(hras2, hierarchies))
    assert(not hras2.refinement_of(hras1, hierarchies))
    assert unicode_strs((hras1, hras2))

def test_refinement_agent_mod_same_order():
    hras1 = Agent('MAP2K1',
        mods=[ModCondition('phosphorylation'),
              ModCondition('phosphorylation')])
    hras2 = Agent('MAP2K1',
        mods=[ModCondition('phosphorylation')])
    assert(hras1.refinement_of(hras2, hierarchies))
    assert(not hras2.refinement_of(hras1, hierarchies))
    assert unicode_strs((hras1, hras2))

def test_refinement_agent_mod_multiple():
    mc1 = ModCondition('phosphorylation', 'S', '218')
    mc2 = ModCondition('phosphorylation', 'S', '298')
    mc3 = ModCondition('phosphorylation', 'S', '222')
    mc4 = ModCondition('phosphorylation')
    mc5 = ModCondition('phosphorylation')

    mek1 = Agent('MAP2K1', mods=[mc1, mc2, mc3])
    mek2 = Agent('MAP2K1', mods=[mc4, mc5])
    erk = Agent('MAPK1')

    st1 = Phosphorylation(mek2, erk)
    st2 = Phosphorylation(mek1, erk, 'T', '185')
    st3 = Phosphorylation(mek1, erk, 'Y', '187')
    assert(st2.refinement_of(st1, hierarchies))
    assert(st3.refinement_of(st1, hierarchies))
    assert(not st1.refinement_of(st2, hierarchies))
    assert(not st1.refinement_of(st3, hierarchies))
    assert unicode_strs((st1, st2, st3))

def test_refinement_agent_mod_generic():
    p = ModCondition('phosphorylation')
    raf3p = Phosphorylation(Agent('RAF', mods=[p,p,p]), Agent('MAP2K1'))
    raf2p = Phosphorylation(Agent('RAF', mods=[p,p]), Agent('MAP2K1'))
    assert(raf3p.refinement_of(raf2p, hierarchies))
    assert(not raf2p.refinement_of(raf3p, hierarchies))
    assert unicode_strs((raf3p, raf2p))

# Check matches implementations for all statement types ---------------------
def test_matches_selfmod():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = Autophosphorylation(nras1, 'tyrosine', '32',
                              evidence=Evidence(text='foo'))
    st2 = Autophosphorylation(nras1, 'tyrosine', '32',
                              evidence=Evidence(text='bar'))
    st3 = Autophosphorylation(nras2, evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))
    assert unicode_strs((st1, st2, st3))

def test_matches_activation():
    """Test matching of entities only, entities match only on name."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = Activation(src, 'kinase1',
                     nras1, 'gtpbound1', True,
                     evidence=Evidence(text='foo'))
    st2 = Activation(src, 'kinase1',
                     nras1, 'gtpbound1', True,
                     evidence=Evidence(text='bar'))
    st3 = Activation(src, 'kinase2',
                     nras2, 'gtpbound2', True,
                     evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))
    assert unicode_strs((st1, st2, st3))

def test_matches_activitymod():
    """Test matching of entities only, entities match only on name."""
    mc = ModCondition('phosphorylation', 'Y', '32')
    mc2 = ModCondition('phosphorylation')
    nras1 = Agent('NRAS', mods=[mc], db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', mods=[mc2], db_refs = {'HGNC': 'dummy'})
    st1 = ActiveForm(nras1, 'gtpbound1', True,
                     evidence=Evidence(text='foo'))
    st2 = ActiveForm(nras1, 'gtpbound1', True,
                     evidence=Evidence(text='bar'))
    st3 = ActiveForm(nras2, 'gtpbound2', True,
                     evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))
    assert unicode_strs((st1, st2, st3))

def test_matches_activatingsub():
    """Test matching of entities only, entities match only on name."""
    mut1 = MutCondition('12', 'G', 'D')
    mut2 = MutCondition('61', 'Q', 'L')
    nras1 = Agent('NRAS', mutations=[mut1], db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', mutations=[mut2], db_refs = {'HGNC': 'dummy'})

    st1 = ActiveForm(nras1, 'gtpbound1', True,
                     evidence=Evidence(text='foo'))
    st2 = ActiveForm(nras1, 'gtpbound1', True,
                     evidence=Evidence(text='bar'))
    st3 = ActiveForm(nras2, 'gtpbound2', True,
                     evidence=Evidence(text='bar'))
    st4 = ActiveForm(nras2, 'gtpbound2', False,
                     evidence=Evidence(text='bar'))
    st5 = ActiveForm(nras2, 'gtpbound3', True,
                     evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))
    assert(not st3.matches(st4)) # Differ only in relationship
    assert(not st3.matches(st5)) # Differ only in activity
    assert unicode_strs((st1, st2, st3, st4, st5))

def test_matches_rasgef():
    """Test matching of entities only, entities match only on name."""
    sos1 = Agent('SOS1', db_refs = {'HGNC': 'sos1'})
    sos2 = Agent('SOS1', db_refs = {'HGNC': 'sos2'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = RasGef(sos1, 'gtpbound1', nras1,
                 evidence=Evidence(text='foo'))
    st2 = RasGef(sos1, 'gtpbound1', nras1,
                 evidence=Evidence(text='bar'))
    st3 = RasGef(sos2, 'gtpbound2', nras2,
                 evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))
    assert unicode_strs((st1, st2, st3))

def test_matches_rasgap():
    rasa1 = Agent('RASA1', db_refs = {'HGNC': 'rasa1'})
    rasa2 = Agent('RASA1', db_refs = {'HGNC': 'rasa2'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = RasGap(rasa1, 'gtpbound1', nras1,
                 evidence=Evidence(text='foo'))
    st2 = RasGap(rasa1, 'gtpbound1', nras1,
                 evidence=Evidence(text='bar'))
    st3 = RasGap(rasa2, 'gtpbound2', nras2,
                 evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))
    assert unicode_strs((st1, st2, st3))

def test_matches_complex():
    ksr1 = Agent('KSR1', db_refs = {'HGNC': 'ksr1'})
    ksr2 = Agent('KSR1', db_refs = {'HGNC': 'ksr2'})
    braf1 = Agent('BRAF', db_refs = {'HGNC': 'braf1'})
    braf2 = Agent('BRAF', db_refs = {'HGNC': 'braf2'})
    map2k1 = Agent('MAP2K1', db_refs = {'HGNC': 'map2k1'})
    map2k2 = Agent('MAP2K1', db_refs = {'HGNC': 'map2k2'})
    st1 = Complex([ksr1, braf1, map2k1], evidence=Evidence(text='foo'))
    st2 = Complex([ksr1, braf1, map2k1], evidence=Evidence(text='bar'))
    st3 = Complex([braf1, map2k1, ksr1], evidence=Evidence(text='bax'))
    assert(st1.matches(st2))
    assert(st2.matches(st3))
    assert(st3.matches(st1))
    assert unicode_strs((st1, st2, st3))

# Entity matching between statements ----------------------------------------
def test_agent_entity_match():
    """Agents match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    assert(nras1.entity_matches(nras2))
    assert unicode_strs((nras1, nras2))

def test_entities_match_mod():
    """Test matching of entities only, entities match only on name."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = Phosphorylation(src, nras1, 'tyrosine', '32',
                          evidence=Evidence(text='foo'))
    st2 = Phosphorylation(src, nras2,
                          evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))
    assert unicode_strs((st1, st2))

def test_entities_match_selfmod():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = Autophosphorylation(nras1, 'tyrosine', '32',
                          evidence=Evidence(text='foo'))
    st2 = Autophosphorylation(nras2,
                          evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))
    assert unicode_strs((st1, st2))

def test_entities_match_activation():
    """Test matching of entities only, entities match only on name."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = Activation(src, 'Kinase1',
                     nras1, 'gtpbound1', True,
                     evidence=Evidence(text='foo'))
    st2 = Activation(src, 'Kinase2',
                     nras2, 'gtpbound2', True,
                     evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))
    assert unicode_strs((st1, st2))

def test_entities_match_activitymod():
    """Test matching of entities only, entities match only on name."""
    mc1 = ModCondition('phosphorylation', 'tyrosine', '32')
    mc2 = ModCondition('phosphorylation')
    nras1 = Agent('NRAS', mods=[mc1], db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', mods=[mc2], db_refs={'HGNC': 'dummy'})
    st1 = ActiveForm(nras1, 'gtpbound1', True,
                     evidence=Evidence(text='foo'))
    st2 = ActiveForm(nras2, 'gtpbound2', False,
                     evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))
    assert unicode_strs((st1, st2))

def test_entities_match_activatingsub():
    """Test matching of entities only, entities match only on name."""
    mc1 = MutCondition('12', 'G', 'D')
    mc2 = MutCondition('61', 'Q', 'L')
    nras1 = Agent('NRAS', mutations=[mc1], db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', mutations=[mc2], db_refs = {'HGNC': 'dummy'})
    st1 = ActiveForm(nras1, 'gtpbound1', True,
                     evidence=Evidence(text='foo'))
    st2 = ActiveForm(nras2, 'gtpbound2', False,
                     evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))
    assert unicode_strs((st1, st2))

def test_entities_match_rasgef():
    """Test matching of entities only, entities match only on name."""
    sos1 = Agent('SOS1', db_refs = {'HGNC': 'sos1'})
    sos2 = Agent('SOS1', db_refs = {'HGNC': 'sos2'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = RasGef(sos1, 'gtpbound1', nras1,
                 evidence=Evidence(text='foo'))
    st2 = RasGef(sos2, 'gtpbound2', nras2,
                 evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))
    assert unicode_strs((st1, st2))

def test_entities_match_rasgap():
    """Test matching of entities only, entities match only on name."""
    rasa1 = Agent('RASA1', db_refs = {'HGNC': 'rasa1'})
    rasa2 = Agent('RASA1', db_refs = {'HGNC': 'rasa2'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = RasGap(rasa1, 'gtpbound1', nras1,
                 evidence=Evidence(text='foo'))
    st2 = RasGap(rasa2, 'gtpbound2', nras2,
                 evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_complex():
    """Test matching of entities only, entities match only on name."""
    ksr1 = Agent('KSR1', db_refs = {'HGNC': 'ksr1'})
    ksr2 = Agent('KSR1', db_refs = {'HGNC': 'ksr2'})
    braf1 = Agent('BRAF', db_refs = {'HGNC': 'braf1'})
    braf2 = Agent('BRAF', db_refs = {'HGNC': 'braf2'})
    map2k1 = Agent('MAP2K1', db_refs = {'HGNC': 'map2k1'})
    map2k2 = Agent('MAP2K1', db_refs = {'HGNC': 'map2k2'})
    st1 = Complex([ksr1, braf1, map2k1], evidence=Evidence(text='foo'))
    st2 = Complex([ksr2, braf2, map2k2], evidence=Evidence(text='bar'))
    st3 = Complex([braf2, map2k2, ksr2], evidence=Evidence(text='baz'))
    assert(st1.entities_match(st2))
    assert(st2.entities_match(st3))
    assert(st3.entities_match(st1))

def test_agent_superfamily_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    ras = Agent('RAS', db_refs = {'BE': 'RAS'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    assert nras.refinement_of(ras, hierarchies)
    assert not ras.refinement_of(nras, hierarchies)
    # The top-level list should contain only one statement, the gene-level
    # one, supported by the family one.

def test_agent_boundcondition_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    bc1 = BoundCondition(Agent('BRAF', db_refs = {'HGNC': '1097'}), True)
    bc2 = BoundCondition(Agent('RAF1', db_refs = {'HGNC': '9829'}), True)
    bc3 = BoundCondition(Agent('RAF1', db_refs = {'HGNC': '9829'}), False)
    bc4 = BoundCondition(Agent('RAF', db_refs = {'BE': 'RAF'}), True)

    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'}, bound_conditions=[bc1])
    nras2 = Agent('NRAS', db_refs = {'HGNC': '7989'}, bound_conditions=[bc2])
    nras3 = Agent('NRAS', db_refs = {'HGNC': '7989'}, bound_conditions=[bc3])
    nras4 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras5 = Agent('NRAS', db_refs = {'HGNC': '7989'},
                  bound_conditions=[bc4])

    # nras1 (bound to BRAF)
    assert not nras2.refinement_of(nras1, hierarchies)
    assert not nras3.refinement_of(nras1, hierarchies)
    assert not nras4.refinement_of(nras1, hierarchies)
    assert not nras5.refinement_of(nras1, hierarchies)
    # nras2 (bound to CRAF)
    assert not nras1.refinement_of(nras2, hierarchies)
    assert not nras3.refinement_of(nras2, hierarchies) # Not bound condition
    assert not nras4.refinement_of(nras2, hierarchies)
    assert not nras5.refinement_of(nras2, hierarchies)
    # nras3 (not bound to CRAF)
    assert not nras1.refinement_of(nras3, hierarchies)
    assert not nras2.refinement_of(nras3, hierarchies) # Not bound condition
    assert not nras4.refinement_of(nras3, hierarchies)
    assert not nras5.refinement_of(nras3, hierarchies)
    # nras4 (no bound condition)
    assert nras1.refinement_of(nras4, hierarchies)
    assert nras2.refinement_of(nras4, hierarchies)
    assert nras3.refinement_of(nras4, hierarchies)
    assert nras5.refinement_of(nras4, hierarchies)
    # nras5 (RAF family bound condition)
    assert nras1.refinement_of(nras5, hierarchies)
    assert nras2.refinement_of(nras5, hierarchies)
    assert not nras3.refinement_of(nras5, hierarchies)
    assert not nras4.refinement_of(nras5, hierarchies)

def test_agent_modification_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    mek1 = Agent('MAP2K1', db_refs = {'HGNC': 'asdf'},
                mods=ModCondition('phosphorylation'))
    mek2 = Agent('MAP2K1', db_refs = {'HGNC': 'asdf'},
                mods=ModCondition('phosphorylation', position='218'))
    mek3 = Agent('MAP2K1', db_refs = {'HGNC': 'asdf'},
                mods=ModCondition('phosphorylation', position='222'))
    mek4 = Agent('MAP2K1', db_refs = {'HGNC': 'asdf'},
                mods=[ModCondition('phosphorylation', position='218'),
                      ModCondition('phosphorylation', position='222')])
    mek5 = Agent('MAP2K1', db_refs = {'HGNC': 'asdf'},
                mods=ModCondition('phosphorylation', 'serine', None))
    mek6 = Agent('MAP2K1', db_refs = {'HGNC': 'asdf'},
                mods=ModCondition('phosphorylation', 'serine', '218'))
    mek7 = Agent('MAP2K1', db_refs = {'HGNC': 'asdf'},
                mods=ModCondition('phosphorylation', 'serine', '222'))
    mek8 = Agent('MAP2K1', db_refs = {'HGNC': 'asdf'},
                mods=[ModCondition('phosphorylation', 'serine', '218'),
                      ModCondition('phosphorylation', 'serine', '222')])

    # mek1 agent is refined by all others
    assert mek2.refinement_of(mek1, hierarchies)
    assert mek3.refinement_of(mek1, hierarchies)
    assert mek4.refinement_of(mek1, hierarchies)
    assert mek5.refinement_of(mek1, hierarchies)
    assert mek6.refinement_of(mek1, hierarchies)
    assert mek7.refinement_of(mek1, hierarchies)
    assert mek8.refinement_of(mek1, hierarchies)
    # mek2
    assert not mek1.refinement_of(mek2, hierarchies)
    assert not mek3.refinement_of(mek2, hierarchies) # Different site
    assert mek4.refinement_of(mek2, hierarchies)
    assert not mek5.refinement_of(mek2, hierarchies) # Cross-relationship
    assert mek6.refinement_of(mek2, hierarchies)
    assert not mek7.refinement_of(mek2, hierarchies) # Different site
    assert mek8.refinement_of(mek2, hierarchies)
    # mek3
    assert not mek1.refinement_of(mek3, hierarchies)
    assert not mek2.refinement_of(mek3, hierarchies)
    assert mek4.refinement_of(mek3, hierarchies)
    assert not mek5.refinement_of(mek3, hierarchies)
    assert not mek6.refinement_of(mek3, hierarchies)
    assert mek7.refinement_of(mek3, hierarchies)
    assert mek8.refinement_of(mek3, hierarchies)
    # mek4
    assert not mek1.refinement_of(mek4, hierarchies)
    assert not mek2.refinement_of(mek4, hierarchies)
    assert not mek3.refinement_of(mek4, hierarchies)
    assert not mek5.refinement_of(mek4, hierarchies)
    assert not mek6.refinement_of(mek4, hierarchies)
    assert not mek7.refinement_of(mek4, hierarchies)
    assert mek8.refinement_of(mek4, hierarchies)
    # mek5
    assert not mek1.refinement_of(mek5, hierarchies)
    assert not mek2.refinement_of(mek5, hierarchies)
    assert not mek3.refinement_of(mek5, hierarchies)
    assert not mek4.refinement_of(mek5, hierarchies)
    assert mek6.refinement_of(mek5, hierarchies)
    assert mek7.refinement_of(mek5, hierarchies)
    assert mek8.refinement_of(mek5, hierarchies)
    # mek6
    assert not mek1.refinement_of(mek6, hierarchies)
    assert not mek2.refinement_of(mek6, hierarchies)
    assert not mek3.refinement_of(mek6, hierarchies)
    assert not mek4.refinement_of(mek6, hierarchies)
    assert not mek5.refinement_of(mek6, hierarchies)
    assert not mek7.refinement_of(mek6, hierarchies)
    assert mek8.refinement_of(mek6, hierarchies)
    # mek7
    assert not mek1.refinement_of(mek7, hierarchies)
    assert not mek2.refinement_of(mek7, hierarchies)
    assert not mek3.refinement_of(mek7, hierarchies)
    assert not mek4.refinement_of(mek7, hierarchies)
    assert not mek5.refinement_of(mek7, hierarchies)
    assert not mek6.refinement_of(mek7, hierarchies)
    assert mek8.refinement_of(mek7, hierarchies)
    # mek8
    assert not mek1.refinement_of(mek8, hierarchies)
    assert not mek2.refinement_of(mek8, hierarchies)
    assert not mek3.refinement_of(mek8, hierarchies)
    assert not mek4.refinement_of(mek8, hierarchies)
    assert not mek5.refinement_of(mek8, hierarchies)
    assert not mek6.refinement_of(mek8, hierarchies)
    assert not mek7.refinement_of(mek8, hierarchies)

def test_phosphorylation_modification_refinement():
    braf = Agent('BRAF', db_refs = {'HGNC': 'braf'})
    mek1 = Agent('MAP2K1', db_refs = {'HGNC': 'map2k1'})
    p1 = Phosphorylation(braf, mek1)
    p2 = Phosphorylation(braf, mek1, position='218')
    p3 = Phosphorylation(braf, mek1, position='222')
    p4 = Phosphorylation(braf, mek1, 'serine')
    p5 = Phosphorylation(braf, mek1, 'serine', '218')
    p6 = Phosphorylation(braf, mek1, 'serine', '222')

    # p1
    assert p2.refinement_of(p1, hierarchies)
    assert p3.refinement_of(p1, hierarchies)
    assert p4.refinement_of(p1, hierarchies)
    assert p5.refinement_of(p1, hierarchies)
    assert p6.refinement_of(p1, hierarchies)
    # p2
    assert not p1.refinement_of(p2, hierarchies)
    assert not p3.refinement_of(p2, hierarchies)
    assert not p4.refinement_of(p2, hierarchies)
    assert p5.refinement_of(p2, hierarchies)
    assert not p6.refinement_of(p2, hierarchies)
    # p3
    assert not p1.refinement_of(p3, hierarchies)
    assert not p2.refinement_of(p3, hierarchies)
    assert not p4.refinement_of(p3, hierarchies)
    assert not p5.refinement_of(p3, hierarchies)
    assert p6.refinement_of(p3, hierarchies)
    # p4
    assert not p1.refinement_of(p4, hierarchies)
    assert not p2.refinement_of(p4, hierarchies)
    assert not p3.refinement_of(p4, hierarchies)
    assert p5.refinement_of(p4, hierarchies)
    assert p6.refinement_of(p4, hierarchies)
    # p5
    assert not p1.refinement_of(p5, hierarchies)
    assert not p2.refinement_of(p5, hierarchies)
    assert not p3.refinement_of(p5, hierarchies)
    assert not p4.refinement_of(p5, hierarchies)
    assert not p6.refinement_of(p5, hierarchies)
    # p6
    assert not p1.refinement_of(p6, hierarchies)
    assert not p2.refinement_of(p6, hierarchies)
    assert not p3.refinement_of(p6, hierarchies)
    assert not p4.refinement_of(p6, hierarchies)
    assert not p5.refinement_of(p6, hierarchies)

def test_autophosphorylation_modification_refinement():
    braf = Agent('BRAF', db_refs = {'HGNC': 'braf'})
    p1 = Autophosphorylation(braf,)
    p2 = Autophosphorylation(braf, position='218')
    p3 = Autophosphorylation(braf, position='222')
    p4 = Autophosphorylation(braf, 'serine')
    p5 = Autophosphorylation(braf, 'serine', '218')
    p6 = Autophosphorylation(braf, 'serine', '222')

    # p1
    assert p2.refinement_of(p1, hierarchies)
    assert p3.refinement_of(p1, hierarchies)
    assert p4.refinement_of(p1, hierarchies)
    assert p5.refinement_of(p1, hierarchies)
    assert p6.refinement_of(p1, hierarchies)
    # p2
    assert not p1.refinement_of(p2, hierarchies)
    assert not p3.refinement_of(p2, hierarchies)
    assert not p4.refinement_of(p2, hierarchies)
    assert p5.refinement_of(p2, hierarchies)
    assert not p6.refinement_of(p2, hierarchies)
    # p3
    assert not p1.refinement_of(p3, hierarchies)
    assert not p2.refinement_of(p3, hierarchies)
    assert not p4.refinement_of(p3, hierarchies)
    assert not p5.refinement_of(p3, hierarchies)
    assert p6.refinement_of(p3, hierarchies)
    # p4
    assert not p1.refinement_of(p4, hierarchies)
    assert not p2.refinement_of(p4, hierarchies)
    assert not p3.refinement_of(p4, hierarchies)
    assert p5.refinement_of(p4, hierarchies)
    assert p6.refinement_of(p4, hierarchies)
    # p5
    assert not p1.refinement_of(p5, hierarchies)
    assert not p2.refinement_of(p5, hierarchies)
    assert not p3.refinement_of(p5, hierarchies)
    assert not p4.refinement_of(p5, hierarchies)
    assert not p6.refinement_of(p5, hierarchies)
    # p6
    assert not p1.refinement_of(p6, hierarchies)
    assert not p2.refinement_of(p6, hierarchies)
    assert not p3.refinement_of(p6, hierarchies)
    assert not p4.refinement_of(p6, hierarchies)
    assert not p5.refinement_of(p6, hierarchies)

def test_activation_modification_refinement():
    raf = Agent('RAF', db_refs={'BE': 'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    mek = Agent('MEK', db_refs={'BE': 'MEK'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})

    st1 = Activation(raf, 'kinase',
                     mek, 'kinase', True)
    st2 = Activation(braf, 'kinase',
                     mek, 'kinase', True)
    st3 = Activation(raf, 'kinase',
                     mek1, 'kinase', True)
    st4 = Activation(braf, 'kinase',
                     mek1, 'kinase', True)
    st5 = Activation(braf, 'kinase',
                     mek1, 'kinase', False)
    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert st4.refinement_of(st1, hierarchies)
    assert not st5.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert st4.refinement_of(st2, hierarchies)
    assert not st5.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st4.refinement_of(st3, hierarchies)
    assert not st5.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)
    assert not st5.refinement_of(st4, hierarchies)
    # st5
    assert not st1.refinement_of(st5, hierarchies)
    assert not st2.refinement_of(st5, hierarchies)
    assert not st3.refinement_of(st5, hierarchies)
    assert not st4.refinement_of(st5, hierarchies)

def test_activation_activity_hierarchy_refinement():
    raf = Agent('RAF', db_refs={'BE': 'RAF'})
    mek = Agent('MEK', db_refs={'BE': 'MEK'})

    st1 = Activation(raf, 'kinase', mek, 'kinase', True)
    st2 = Activation(raf, 'kinase', mek, 'kinase', False)
    st3 = Activation(raf, 'catalytic', mek, 'kinase', True)
    st4 = Activation(raf, 'kinase', mek, 'catalytic', True)
    st5 = Activation(raf, 'catalytic', mek, 'activity', True)
    st6 = Activation(raf, 'activity', mek, 'activity', True)

    assert(not st1.refinement_of(st2, hierarchies))
    assert(not st2.refinement_of(st1, hierarchies))
    assert(st1.refinement_of(st3, hierarchies))
    assert(st1.refinement_of(st4, hierarchies))
    assert(st5.refinement_of(st6, hierarchies))
    assert(st1.refinement_of(st6, hierarchies))
    assert(not st3.refinement_of(st4, hierarchies))
    assert(not st4.refinement_of(st3, hierarchies))

def test_activitymod_refinement():
    mc1 = ModCondition('phosphorylation')
    mc2 = ModCondition('phosphorylation', 'S')
    mc3 = ModCondition('phosphorylation', 'S', '218')
    mc4 = ModCondition('phosphorylation', 'S', '222')
    mek_fam = Agent('MEK')
    mek1 = Agent('MAP2K1')
    p1 = ActiveForm(Agent('MEK', mods=[mc1], db_refs={'BE':'MEK'}),
                    'kinase', True)
    p2 = ActiveForm(Agent('MEK', mods=[mc3], db_refs={'BE':'MEK'}),
                    'kinase', True)
    p3 = ActiveForm(Agent('MAP2K1', mods=[mc1], db_refs={'HGNC':'6840'}),
                    'kinase', True)
    p4 = ActiveForm(Agent('MAP2K1', mods=[mc2], db_refs={'HGNC':'6840'}),
                    'kinase', True)
    p5 = ActiveForm(Agent('MAP2K1', mods=[mc3], db_refs={'HGNC':'6840'}),
                    'kinase', True)
    p6 = ActiveForm(Agent('MAP2K1', mods=[mc4], db_refs={'HGNC':'6840'}),
                    'kinase', True)
    p7 = ActiveForm(Agent('MAP2K1', mods=[mc3, mc4], db_refs={'HGNC':'6840'}),
                    'kinase', True)
    # p1
    assert p2.refinement_of(p1, hierarchies)
    assert p3.refinement_of(p1, hierarchies)
    assert p4.refinement_of(p1, hierarchies)
    assert p5.refinement_of(p1, hierarchies)
    assert p6.refinement_of(p1, hierarchies)
    assert p7.refinement_of(p1, hierarchies)
    # p2
    assert not p1.refinement_of(p2, hierarchies)
    assert not p3.refinement_of(p2, hierarchies)
    assert not p4.refinement_of(p2, hierarchies)
    assert p5.refinement_of(p2, hierarchies)
    assert not p6.refinement_of(p2, hierarchies)
    assert p7.refinement_of(p2, hierarchies)
    # p3
    assert not p1.refinement_of(p3, hierarchies)
    assert not p2.refinement_of(p3, hierarchies)
    assert p4.refinement_of(p3, hierarchies)
    assert p5.refinement_of(p3, hierarchies)
    assert p6.refinement_of(p3, hierarchies)
    assert p7.refinement_of(p3, hierarchies)
    # p4
    assert not p1.refinement_of(p4, hierarchies)
    assert not p2.refinement_of(p4, hierarchies)
    assert not p3.refinement_of(p4, hierarchies)
    assert p5.refinement_of(p4, hierarchies)
    assert p6.refinement_of(p4, hierarchies)
    assert p7.refinement_of(p4, hierarchies)
    # p5
    assert not p1.refinement_of(p5, hierarchies)
    assert not p2.refinement_of(p5, hierarchies)
    assert not p3.refinement_of(p5, hierarchies)
    assert not p4.refinement_of(p5, hierarchies)
    assert not p6.refinement_of(p5, hierarchies)
    assert p7.refinement_of(p5, hierarchies)
    # p6
    assert not p1.refinement_of(p6, hierarchies)
    assert not p2.refinement_of(p6, hierarchies)
    assert not p3.refinement_of(p6, hierarchies)
    assert not p4.refinement_of(p6, hierarchies)
    assert not p5.refinement_of(p6, hierarchies)
    assert p7.refinement_of(p6, hierarchies)
    # p7
    assert not p1.refinement_of(p7, hierarchies)
    assert not p2.refinement_of(p7, hierarchies)
    assert not p3.refinement_of(p7, hierarchies)
    assert not p4.refinement_of(p7, hierarchies)
    assert not p5.refinement_of(p7, hierarchies)
    assert not p6.refinement_of(p7, hierarchies)

def test_activeform_activity_hierarchy_refinement():
    p1 = ActiveForm(Agent('MEK'), 'kinase', True)
    p2 = ActiveForm(Agent('MEK'), 'kinase', False)
    p3 = ActiveForm(Agent('MEK'), 'catalytic', True)
    p4 = ActiveForm(Agent('MEK'), 'activity', True)

    assert(not p1.refinement_of(p2, hierarchies))
    assert(p1.refinement_of(p3, hierarchies))
    assert(p1.refinement_of(p4, hierarchies))
    assert(p3.refinement_of(p4, hierarchies))
    assert(not p4.refinement_of(p3, hierarchies))

def test_activatingsub_family_refinement():
    mc = MutCondition('12', 'G', 'D')
    ras = Agent('RAS', mutations=[mc], db_refs={'BE':'RAS'})
    kras = Agent('KRAS', mutations=[mc], db_refs={'HGNC':'6407'})
    nras = Agent('NRAS', mutations=[mc], db_refs={'HGNC':'7989'})
    st1 = ActiveForm(ras, 'activity', True)
    st2 = ActiveForm(kras, 'activity', True)
    st3 = ActiveForm(nras, 'activity', True)
    st4 = ActiveForm(kras, 'activity', False)
    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert not st4.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert not st4.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert not st4.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)

def test_rasgef_family_refinement():
    sos = Agent('SOS', db_refs={'BE':'SOS'})
    sos1 = Agent('SOS1', db_refs={'HGNC':'11187'})
    ras = Agent('RAS', db_refs={'BE':'RAS'})
    kras = Agent('KRAS', db_refs={'HGNC':'6407'})
    # Statements
    st1 = RasGef(sos, 'activity', ras)
    st2 = RasGef(sos1, 'activity', ras)
    st3 = RasGef(sos, 'activity', kras)
    st4 = RasGef(sos1, 'activity', kras)
    st5 = RasGef(sos1, 'different_activity', kras)
    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert st4.refinement_of(st1, hierarchies)
    assert not st5.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert st4.refinement_of(st2, hierarchies)
    assert not st5.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st4.refinement_of(st3, hierarchies)
    assert not st5.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)
    assert not st5.refinement_of(st4, hierarchies)
    # st5
    assert not st1.refinement_of(st5, hierarchies)
    assert not st2.refinement_of(st5, hierarchies)
    assert not st3.refinement_of(st5, hierarchies)
    assert not st4.refinement_of(st5, hierarchies)

def test_rasgap_family_refinement():
    rasa = Agent('RASA', db_refs={'BE':'RASA'})
    rasa1 = Agent('RASA1', db_refs={'HGNC':'9871'})
    ras = Agent('RAS', db_refs={'BE':'RAS'})
    kras = Agent('KRAS', db_refs={'HGNC':'6407'})
    # Statements
    st1 = RasGap(rasa, 'activity', ras)
    st2 = RasGap(rasa1, 'activity', ras)
    st3 = RasGap(rasa, 'activity', kras)
    st4 = RasGap(rasa1, 'activity', kras)
    st5 = RasGap(rasa1, 'different_activity', kras)
    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert st4.refinement_of(st1, hierarchies)
    assert not st5.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert st4.refinement_of(st2, hierarchies)
    assert not st5.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st4.refinement_of(st3, hierarchies)
    assert not st5.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)
    assert not st5.refinement_of(st4, hierarchies)
    # st5
    assert not st1.refinement_of(st5, hierarchies)
    assert not st2.refinement_of(st5, hierarchies)
    assert not st3.refinement_of(st5, hierarchies)
    assert not st4.refinement_of(st5, hierarchies)

def test_complex_family_refinement():
    raf = Agent('RAF', db_refs={'BE':'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC':'1097'})
    raf1 = Agent('RAF1', db_refs={'HGNC':'9829'})
    mek = Agent('MEK', db_refs={'BE':'MEK'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC':'6840'})

    st1 = Complex([raf, mek])
    st2 = Complex([braf, mek])
    st3 = Complex([mek1, raf])
    st4 = Complex([braf, mek1])
    st5 = Complex([braf, raf1])

    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert st4.refinement_of(st1, hierarchies)
    assert not st5.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert st4.refinement_of(st2, hierarchies)
    assert not st5.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st4.refinement_of(st3, hierarchies)
    assert not st5.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)
    assert not st5.refinement_of(st4, hierarchies)
    # st5
    assert not st1.refinement_of(st5, hierarchies)
    assert not st2.refinement_of(st5, hierarchies)
    assert not st3.refinement_of(st5, hierarchies)
    assert not st4.refinement_of(st5, hierarchies)

def test_mismatched_complex_refinement():
    ras = Agent('RAS')
    raf = Agent('RAF')
    mek = Agent('MEK')
    st1 = Complex([ras, raf])
    st2 = Complex([mek, ras, raf])
    assert not st1.refinement_of(st2, hierarchies)
    assert not st2.refinement_of(st1, hierarchies)

@raises(InvalidResidueError)
def test_residue_mod_condition():
    mc = ModCondition('phosphorylation', 'xyz')

@raises(InvalidResidueError)
def test_residue_mod():
    Phosphorylation(Agent('a'), Agent('b'), 'xyz')

@raises(InvalidResidueError)
def test_residue_selfmod():
    Autophosphorylation(Agent('a'), 'xyz')

def test_valid_mod_residue():
    mc = ModCondition('phosphorylation', 'serine')
    assert(mc.residue == 'S')
    assert unicode_strs(mc)

def test_valid_residue():
    assert(get_valid_residue('serine') == 'S')
    assert(get_valid_residue('ser') == 'S')
    assert(get_valid_residue('Serine') == 'S')
    assert(get_valid_residue('SERINE') == 'S')

def test_modcondition_order_actmod():
    mc1 = ModCondition('phoshporylation', 'S', '222')
    mc2 = ModCondition('phoshporylation', 'S', '224')
    p1 = ActiveForm(Agent('MAP2K1', mods=[mc1, mc2]), 'kinase', True)
    p2 = ActiveForm(Agent('MAP2K1', mods=[mc2, mc1]), 'kinase', True)
    assert(p1.matches(p2))
    assert unicode_strs((p1, p2))

def test_modcondition_order_agent():
    mc1 = ModCondition('phoshporylation', 'S', '222')
    mc2 = ModCondition('phoshporylation', 'S', '224')
    p1 = Agent('MAP2K1', mods=[mc1, mc2])
    p2 = Agent('MAP2K1', mods=[mc2, mc1])
    assert(p1.matches(p2))
    assert unicode_strs((p1, p2))

def test_eq_mut():
    assert(MutCondition('600', 'V', 'E').equals(MutCondition('600', 'V', 'E')))
    assert(not MutCondition('600', 'V', 'E').equals(
                                             MutCondition('600', 'V', 'D')))

def test_eq_agent():
    assert(Agent('one').equals(Agent('one')))
    assert(not Agent('one').equals(Agent('two')))
    assert(not Agent('one', db_refs={'UP': '123'}).equals(
           Agent('one', db_refs={'UP': '999'})))
    assert(Agent('one', mods=[ModCondition('phosphorylation')]).equals(
           Agent('one', mods=[ModCondition('phosphorylation')])))
    assert(not Agent('one', mods=[ModCondition('phosphorylation')]).equals(
           Agent('one', mods=[ModCondition('ubiquitination')])))
    assert(Agent('one', mutations=[MutCondition('600', 'V', 'E')]).equals(
           Agent('one', mutations=[MutCondition('600', 'V', 'E')])))
    assert(not Agent('one', mutations=[MutCondition('600', 'V', 'E')]).equals(
           Agent('one', mutations=[MutCondition('600', 'V', 'D')])))
    assert(Agent('one',
                 bound_conditions=[BoundCondition(Agent('two'), True)]).equals(
           Agent('one',
                 bound_conditions=[BoundCondition(Agent('two'), True)])))
    assert(not Agent('one',
                     bound_conditions=[BoundCondition(Agent('two'),
                                                      True)]).equals(
           Agent('one', bound_conditions=[BoundCondition(Agent('two'),
                                                         False)])))
    assert(not Agent('one', bound_conditions=[BoundCondition(Agent('two'),
                                                             True)]).equals(
           Agent('one', bound_conditions=[BoundCondition(Agent('three'),
                                                         True)])))

def test_eq_stmt():
    ev1 = Evidence(text='1')
    ev2 = Evidence(text='2')
    assert(Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
            Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1])))
    assert(not Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
            Phosphorylation(Agent('a'), Agent('b'), evidence=[ev2])))
    assert(not Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
            Phosphorylation(Agent('a'), Agent('c'), evidence=[ev2])))
    assert(not Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
            Phosphorylation(Agent('a'), Agent('b'), 'S', evidence=[ev2])))
    assert(Complex([Agent('a'), Agent('b')], evidence=[ev1]).equals(
           Complex([Agent('a'), Agent('b')], evidence=[ev1])))
    assert(not Complex([Agent('a'), Agent('b')], evidence=[ev1]).equals(
           Complex([Agent('a'), Agent('b')], evidence=[ev2])))
    assert(Activation(Agent('a'), 'activity',
                      Agent('b'), 'activity', True, evidence=[ev1]).equals(
           Activation(Agent('a'), 'activity',
                      Agent('b'), 'activity', True, evidence=[ev1])))
    assert(not Activation(Agent('a'), 'activity',
                          Agent('b'), 'activity', True, evidence=[ev1]).equals(
           Activation(Agent('a'), 'activity',
                      Agent('c'), 'activity', True, evidence=[ev1])))
    assert(not Activation(Agent('a'), 'activity',
                          Agent('b'), 'activity', True, evidence=[ev1]).equals(
           Activation(Agent('a'), 'activity',
                      Agent('b'), 'kinase', True, evidence=[ev1])))
    assert(not Activation(Agent('a'), 'activity',
                          Agent('b'), 'activity', True, evidence=[ev1]).equals(
           Activation(Agent('a'), 'activity',
                      Agent('b'), 'activity', True, evidence=[ev2])))

def test_serialize():
    ev1 = Evidence(text='1\U0001F4A9')
    st = Phosphorylation(Agent('a\U0001F4A9'), Agent('b'), evidence=[ev1])
    jstr = st.to_json()
    st2 = Phosphorylation.from_json(jstr)
    assert(st.equals(st2))
    assert unicode_strs((ev1, st, st2))

def test_serialize_errors():
    st = Phosphorylation(Agent('a\U0001F4A9'), Agent('b\U0001F4A9'))
    jstr = st.to_json()
    st2 = Complex.from_json(jstr)
    assert(st2 is None)
    st3 = Phosphorylation.from_json('{}')
    assert(st3 is None)
    st4 = Phosphorylation.from_json('xyz' + jstr)
    assert(st4 is None)
    assert unicode_strs((st, st2, st3, st4))

def test_location_refinement():
    a1 = Agent('a', location='plasma membrane')
    a2 = Agent('a', location='cell')
    a3 = Agent('a', location='cytoplasm')
    a4 = Agent('a')
    a5 = Agent('a')

    assert(a1.refinement_of(a2, hierarchies))
    assert(not a2.refinement_of(a3, hierarchies))
    assert(a4.refinement_of(a5, hierarchies))
    assert(not a1.refinement_of(a3, hierarchies))
    assert(not a3.refinement_of(a1, hierarchies))
    assert(a2.refinement_of(a4, hierarchies))
    assert(a3.refinement_of(a4, hierarchies))

def test_activity_refinement():
    a1 = Agent('a', active='kinase')
    a2 = Agent('a', active='activity')
    a3 = Agent('a', active='catalytic')
    a4 = Agent('a')

    assert(a1.refinement_of(a2, hierarchies))
    assert(not a2.refinement_of(a3, hierarchies))
    assert(not a4.refinement_of(a1, hierarchies))
    assert(a1.refinement_of(a3, hierarchies))
    assert(a3.refinement_of(a2, hierarchies))
    assert(not a3.refinement_of(a1, hierarchies))
    assert(a1.refinement_of(a4, hierarchies))
    assert(a2.refinement_of(a4, hierarchies))

def test_translocation_refinement():
    st1 = Translocation(Agent('a'), 'plasma membrane', 'cytoplasm')
    st2 = Translocation(Agent('a'), 'plasma membrane', None)
    st3 = Translocation(Agent('a'), None, 'cytoplasm')
    st4 = Translocation(Agent('a'), 'cell', 'cytoplasm')
    st5 = Translocation(Agent('a'), 'cell', 'cell')
    st6 = Translocation(Agent('a'), 'plasma membrane', 'cell')
    st7 = Translocation(Agent('a'), 'nucleus', 'cytoplasm')
    st8 = Translocation(Agent('a'), None, 'cell')
    st9 = Translocation(Agent('a'), None, None)
    assert(st3.refinement_of(st8, hierarchies))
    assert(st1.refinement_of(st2, hierarchies))
    assert(st1.refinement_of(st3, hierarchies))
    assert(not st2.refinement_of(st3, hierarchies))
    assert(st1.refinement_of(st4, hierarchies))
    assert(not st2.refinement_of(st4, hierarchies))
    assert(st4.refinement_of(st5, hierarchies))
    assert(st6.refinement_of(st5, hierarchies))
    assert(not st1.refinement_of(st7, hierarchies))
    assert(st7.refinement_of(st4, hierarchies))
    assert(st8.refinement_of(st9, hierarchies))
    assert(st7.refinement_of(st9, hierarchies))

def test_complex_refinement_order():
    st1 = Complex([Agent('MED23'), Agent('ELK1')])
    st2 = Complex([Agent('ELK1', mods=[ModCondition('phosphorylation')]),
                   Agent('MED23')])
    assert(st2.refinement_of(st1, hierarchies))
    assert(not st1.refinement_of(st2, hierarchies))

def test_homodimer_bound_to():
    KRAS = Agent('KRAS')
    HRAS = Agent('HRAS')
    NRAS = Agent('NRAS')
    BRAFK = Agent('BRAF', bound_conditions=[BoundCondition(KRAS, True)])
    BRAFH = Agent('BRAF', bound_conditions=[BoundCondition(HRAS, True)])
    BRAFN = Agent('BRAF', bound_conditions=[BoundCondition(NRAS, True)])

    st1 = Complex([BRAFK, BRAFN])
    st2 = Complex([BRAFN, BRAFK])
    st3 = Complex([BRAFK, BRAFH])
    assert st1.matches(st2)
    assert st2.matches(st1)
    assert not st1.matches(st3)
    assert not st3.matches(st2)

def test_unicode_str_methods():
    ag = Agent('MAPK1\U0001F4A9')
    print(ag)
    ev = Evidence(text='foo \U0001F4A9 bar')
    print(ev)
    print(repr(ev))

    st = Phosphorylation(ag, ag, evidence=ev)
    print(st)
    print(repr(st))

    st1 = Autophosphorylation(ag, evidence=ev)
    print(st1)
    print(repr(st1))

    st = Activation(ag, 'activity', ag, 'activity', is_activation=True,
                    evidence=ev)
    print(st)
    print(repr(st))

    st = ActiveForm(ag, 'activity', True)
    print(st)
    print(repr(st))

    st = HasActivity(ag, 'activity', True)
    print(st)
    print(repr(st))

    st = RasGef(ag, 'gef', ag, evidence=ev)
    print(st)
    print(repr(st))

    st = RasGap(ag, 'gap', ag, evidence=ev)
    print(st)
    print(repr(st))

    st = Complex([ag, ag], evidence=ev)
    print(st)
    print(repr(st))

