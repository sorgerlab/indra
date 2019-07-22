from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import json
import unittest
from copy import deepcopy
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
    assert st1.matches(st2)
    assert unicode_strs(st1)


def test_matches_key():
    ras = Agent('Ras')
    raf = Agent('Raf')
    st1 = Complex([ras, raf])
    st2 = Complex([ras, raf])
    assert st1.matches_key() == st2.matches_key()
    assert unicode_strs(st1)


def test_matches_key_unicode():
    ras = Agent('Ras')
    rasu = Agent(u'Ras')
    raf = Agent('Raf')
    st1 = Complex([ras, raf])
    st2 = Complex([rasu, raf])
    assert st1.matches_key() == st2.matches_key()
    assert unicode_strs(st1)
    assert unicode_strs(st2)


def test_matches_key_unicode2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek, u'S')
    st2 = Phosphorylation(raf, mek, 'S')
    assert st1.matches_key() == st2.matches_key()
    assert unicode_strs(st1)
    assert unicode_strs(st2)


def test_matches_key_unicode3():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek, 'S', u'222')
    st2 = Phosphorylation(raf, mek, 'S', '222')
    assert st1.matches_key() == st2.matches_key()
    assert unicode_strs(st1)
    assert unicode_strs(st2)


def test_matches2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek)
    assert st1.matches(st2)
    assert unicode_strs(st1)
    assert unicode_strs(st2)


def test_matches_key2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek)
    assert st1.matches_key() == st2.matches_key()
    assert unicode_strs(st1)
    assert unicode_strs(st2)


def test_not_matches():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek, 'tyrosine')
    assert not st1.matches(st2)
    assert unicode_strs(st1)
    assert unicode_strs(st2)


def test_not_matches_key():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek, 'tyrosine')
    assert st1.matches_key() != st2.matches_key()
    assert unicode_strs(st1)
    assert unicode_strs(st2)


def test_matches_dbrefs():
    hras1 = Agent('HRAS', db_refs={'hgnc': 1111})
    hras2 = Agent('HRAS', db_refs={'hgnc': 9999})
    assert hras1.matches(hras2)
    assert unicode_strs(hras1)
    assert unicode_strs(hras2)


def test_matches_key_dbrefs():
    hras1 = Agent('HRAS', db_refs={'hgnc': 1111})
    hras2 = Agent('HRAS', db_refs={'hgnc': 9999})
    assert hras1.matches_key() == hras2.matches_key()
    assert unicode_strs((hras1, hras2))


def test_matches_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    assert hras1.matches(hras2)
    assert unicode_strs((hras1, hras2))


def test_matches_key_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    assert hras1.matches_key() == hras2.matches_key()
    assert unicode_strs((hras1, hras2))


def test_not_matches_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('RAF1'), True)])
    assert not hras1.matches(hras2)
    assert unicode_strs((hras1, hras2))


def test_not_matches_key_bound():
    hras1 = Agent('HRAS',
                  bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
                  bound_conditions=[BoundCondition(Agent('RAF1'), True)])
    assert hras1.matches_key() != hras2.matches_key()
    assert unicode_strs((hras1, hras2))


def test_not_matches_bound2():
    hras1 = Agent('HRAS',
                  bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
                  bound_conditions=[BoundCondition(Agent('BRAF'), False)])
    assert not hras1.matches(hras2)
    assert unicode_strs((hras1, hras2))


def test_not_matches_key_bound2():
    hras1 = Agent('HRAS',
                  bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
                  bound_conditions=[BoundCondition(Agent('BRAF'), False)])
    assert hras1.matches_key() != hras2.matches_key()
    assert unicode_strs((hras1, hras2))


def test_matches_bound_multiple():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    assert hras1.matches(hras2)
    assert unicode_strs((hras1, hras2))


def test_matches_key_bound_multiple():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    assert hras1.matches_key() == hras2.matches_key()
    assert unicode_strs((hras1, hras2))


def test_matches_bound_multiple_order():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('RAF1'), True),
                                            BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    assert hras1.matches(hras2)
    assert unicode_strs((hras1, hras2))


def test_matches_key_bound_multiple_order():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('RAF1'), True),
                                            BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    assert hras1.matches_key() == hras2.matches_key()
    assert unicode_strs((hras1, hras2))


def test_matches_agent_mod_order():
    hras1 = Agent('MAP2K1',
        mods=[ModCondition('phosphorylation'), ModCondition('ubiquitination')])
    hras2 = Agent('MAP2K1',
        mods=[ModCondition('ubiquitination'), ModCondition('phosphorylation')])
    assert hras1.matches(hras2)
    assert unicode_strs((hras1, hras2))


def test_refinement_agent_mod_order():
    hras1 = Agent('MAP2K1',
                  mods=[ModCondition('phosphorylation', 'S'),
                        ModCondition('ubiquitination')])
    hras2 = Agent('MAP2K1',
                  mods=[ModCondition('ubiquitination'),
                        ModCondition('phosphorylation')])
    assert hras1.refinement_of(hras2, hierarchies)
    assert not hras2.refinement_of(hras1, hierarchies)
    assert unicode_strs((hras1, hras2))


def test_refinement_agent_mod_same_order():
    hras1 = Agent('MAP2K1',
                  mods=[ModCondition('phosphorylation'),
                        ModCondition('phosphorylation')])
    hras2 = Agent('MAP2K1',
                  mods=[ModCondition('phosphorylation')])
    assert hras1.refinement_of(hras2, hierarchies)
    assert not hras2.refinement_of(hras1, hierarchies)
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
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert not st1.refinement_of(st2, hierarchies)
    assert not st1.refinement_of(st3, hierarchies)
    assert unicode_strs((st1, st2, st3))


def test_refinement_agent_mod_generic():
    p = ModCondition('phosphorylation')
    raf3p = Phosphorylation(Agent('RAF', mods=[p,p,p]), Agent('MAP2K1'))
    raf2p = Phosphorylation(Agent('RAF', mods=[p,p]), Agent('MAP2K1'))
    assert raf3p.refinement_of(raf2p, hierarchies)
    assert not raf2p.refinement_of(raf3p, hierarchies)
    assert unicode_strs((raf3p, raf2p))


# Check matches implementations for all statement types ---------------------
def test_matches_selfmod():
    """Test matching of entities only."""
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    st1 = Autophosphorylation(nras1, 'tyrosine', '32',
                              evidence=Evidence(text='foo'))
    st2 = Autophosphorylation(nras1, 'tyrosine', '32',
                              evidence=Evidence(text='bar'))
    st3 = Autophosphorylation(nras2, evidence=Evidence(text='bar'))
    assert st1.matches(st2)
    assert not st1.matches(st3)
    assert unicode_strs((st1, st2, st3))


def test_matches_activation():
    """Test matching of entities only."""
    src = Agent('SRC', db_refs={'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    st1 = Activation(src, nras1, 'gtpbound',
                     evidence=Evidence(text='foo'))
    st2 = Activation(src, nras1, 'gtpbound',
                     evidence=Evidence(text='bar'))
    st3 = Activation(src, nras2, 'phosphatase',
                     evidence=Evidence(text='bar'))
    st4 = Inhibition(src, nras2, 'phosphatase',
                     evidence=Evidence(text='bar'))
    assert st1.matches(st2)
    assert not st1.matches(st3)
    assert not st3.matches(st4)
    assert unicode_strs((st1, st2, st3))


def test_matches_activitymod():
    """Test matching of entities only."""
    mc = ModCondition('phosphorylation', 'Y', '32')
    mc2 = ModCondition('phosphorylation')
    nras1 = Agent('NRAS', mods=[mc], db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', mods=[mc2], db_refs={'HGNC': 'dummy'})
    st1 = ActiveForm(nras1, 'gtpbound', True,
                     evidence=Evidence(text='foo'))
    st2 = ActiveForm(nras1, 'gtpbound', True,
                     evidence=Evidence(text='bar'))
    st3 = ActiveForm(nras2, 'phosphatase', True,
                     evidence=Evidence(text='bar'))
    assert st1.matches(st2)
    assert not st1.matches(st3)
    assert unicode_strs((st1, st2, st3))


def test_matches_activatingsub():
    """Test matching of entities only."""
    mut1 = MutCondition('12', 'G', 'D')
    mut2 = MutCondition('61', 'Q', 'L')
    nras1 = Agent('NRAS', mutations=[mut1], db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', mutations=[mut2], db_refs={'HGNC': 'dummy'})

    st1 = ActiveForm(nras1, 'gtpbound', True,
                     evidence=Evidence(text='foo'))
    st2 = ActiveForm(nras1, 'gtpbound', True,
                     evidence=Evidence(text='bar'))
    st3 = ActiveForm(nras2, 'phosphatase', True,
                     evidence=Evidence(text='bar'))
    st4 = ActiveForm(nras2, 'phosphatase', False,
                     evidence=Evidence(text='bar'))
    st5 = ActiveForm(nras2, 'kinase', True,
                     evidence=Evidence(text='bar'))
    assert st1.matches(st2)
    assert not st1.matches(st3)
    assert not st3.matches(st4) # Differ only in relationship
    assert not st3.matches(st5) # Differ only in activity
    assert unicode_strs((st1, st2, st3, st4, st5))


def test_matches_gef():
    """Test matching of entities only."""
    sos1 = Agent('SOS1', db_refs={'HGNC': 'sos1'})
    sos2 = Agent('SOS1', db_refs={'HGNC': 'sos2'})
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    st1 = Gef(sos1, nras1,
              evidence=Evidence(text='foo'))
    st2 = Gef(sos1, nras1,
              evidence=Evidence(text='bar'))
    st3 = Gef(sos2, nras2,
              evidence=Evidence(text='bar'))
    assert st1.matches(st2)
    assert not st1.matches(st3)
    assert unicode_strs((st1, st2, st3))


def test_matches_gap():
    rasa1 = Agent('RASA1', db_refs={'HGNC': 'rasa1'})
    rasa2 = Agent('RASA1', db_refs={'HGNC': 'rasa2'})
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    st1 = Gap(rasa1, nras1,
              evidence=Evidence(text='foo'))
    st2 = Gap(rasa1, nras1,
              evidence=Evidence(text='bar'))
    st3 = Gap(rasa2, nras2,
              evidence=Evidence(text='bar'))
    assert st1.matches(st2)
    assert not st1.matches(st3)
    assert unicode_strs((st1, st2, st3))


def test_matches_complex():
    ksr1 = Agent('KSR1', db_refs={'HGNC': 'ksr1'})
    ksr2 = Agent('KSR1', db_refs={'HGNC': 'ksr2'})
    braf1 = Agent('BRAF', db_refs={'HGNC': 'braf1'})
    braf2 = Agent('BRAF', db_refs={'HGNC': 'braf2'})
    map2k1 = Agent('MAP2K1', db_refs={'HGNC': 'map2k1'})
    map2k2 = Agent('MAP2K1', db_refs={'HGNC': 'map2k2'})
    st1 = Complex([ksr1, braf1, map2k1], evidence=Evidence(text='foo'))
    st2 = Complex([ksr1, braf1, map2k1], evidence=Evidence(text='bar'))
    st3 = Complex([braf1, map2k1, ksr1], evidence=Evidence(text='bax'))
    assert st1.matches(st2)
    assert st2.matches(st3)
    assert st3.matches(st1)
    assert unicode_strs((st1, st2, st3))


# Entity matching between statements ----------------------------------------
def test_agent_entity_match():
    """Agents match on name and grounding."""
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras3 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    assert nras1.entity_matches(nras2)
    assert not nras1.entity_matches(nras3)
    assert unicode_strs((nras1, nras2, nras3))


def test_agent_entity_match_chebi():
    vem1 = Agent('vemurafenib', db_refs={'CHEBI': 'CHEBI:63637'})
    vem2 = Agent('Vemurafenib', db_refs={'CHEBI': 'CHEBI:63637',
                                         'PUBCHEM': 'XXX'})
    vem3 = Agent('vemurafenib', db_refs={'CHEBI': 'XXX',
                                         'PUBCHEM': 'XXX'})
    assert vem1.entity_matches(vem2)
    assert not vem1.entity_matches(vem3)
    assert not vem2.entity_matches(vem3)


def test_agent_entity_match_go_mesh():
    adh1 = Agent('ADHESION', db_refs={'GO': 'GO:0007155',
                                      'MESH': 'D002448'})
    adh2 = Agent('adhesion', db_refs={'GO': 'GO:0007155'})
    adh3 = Agent('Adhesion', db_refs={'MESH': 'D002448'})
    # These are satisfied because GO takes priority over MESH
    assert adh1.entity_matches(adh2)
    assert not adh1.entity_matches(adh3)
    assert not adh2.entity_matches(adh3)


def test_agent_entity_match_ungrounded():
    ag1 = Agent('something', db_refs={'XXX': 'XXX'})
    ag2 = Agent('something', db_refs={'YYY': 'YYY'})
    ag3 = Agent('Something', db_refs={'YYY': 'YYY'})
    ag4 = Agent('SOMETHING', db_refs={'ZZZ': 'ZZZ'})
    assert ag1.entity_matches(ag2)
    assert not ag1.entity_matches(ag3)
    assert not ag1.entity_matches(ag4)
    assert not ag2.entity_matches(ag3)
    assert not ag3.entity_matches(ag4)


def test_entities_match_mod():
    """Test matching of entities only, entities match on name and grounding."""
    src = Agent('SRC', db_refs={'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras3 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    st1 = Phosphorylation(src, nras1, 'tyrosine', '32',
                          evidence=Evidence(text='foo'))
    st2 = Phosphorylation(src, nras2,
                          evidence=Evidence(text='bar'))
    st3 = Phosphorylation(src, nras3,
                          evidence=Evidence(text='baz'))
    assert st1.entities_match(st2)
    assert not st1.entities_match(st3)
    assert unicode_strs((st1, st2, st3))


def test_entities_match_selfmod():
    """Test matching of entities only, entities match on name and grounding."""
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras3 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    st1 = Autophosphorylation(nras1, 'tyrosine', '32',
                              evidence=Evidence(text='foo'))
    st2 = Autophosphorylation(nras2,
                              evidence=Evidence(text='bar'))
    st3 = Autophosphorylation(nras3,
                              evidence=Evidence(text='baz'))
    assert st1.entities_match(st2)
    assert not st1.entities_match(st3)
    assert unicode_strs((st1, st2, st3))


def test_entities_match_activation():
    """Test matching of entities only, entities match on name and grounding."""
    src = Agent('SRC', db_refs={'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras3 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    st1 = Activation(src, nras1, 'gtpbound',
                     evidence=Evidence(text='foo'))
    st2 = Activation(src, nras2, 'phosphatase',
                     evidence=Evidence(text='bar'))
    st3 = Activation(src, nras3, 'phosphatase',
                     evidence=Evidence(text='baz'))
    assert st1.entities_match(st2)
    assert not st1.entities_match(st3)
    assert unicode_strs((st1, st2, st3))


def test_entities_match_activitymod():
    """Test matching of entities only, entities match on name and grounding."""
    mc1 = ModCondition('phosphorylation', 'tyrosine', '32')
    mc2 = ModCondition('phosphorylation')
    nras1 = Agent('NRAS', mods=[mc1], db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', mods=[mc2], db_refs={'HGNC': '7989'})
    nras3 = Agent('NRAS', mods=[mc1], db_refs={'HGNC': 'dummy'})
    st1 = ActiveForm(nras1, 'gtpbound', True,
                     evidence=Evidence(text='foo'))
    st2 = ActiveForm(nras2, 'phosphatase', False,
                     evidence=Evidence(text='bar'))
    st3 = ActiveForm(nras3, 'gtpbound', False,
                     evidence=Evidence(text='baz'))
    assert st1.entities_match(st2)
    assert not st1.entities_match(st3)
    assert unicode_strs((st1, st2, st3))


def test_entities_match_activatingsub():
    """Test matching of entities only, entities match on name and grounding."""
    mc1 = MutCondition('12', 'G', 'D')
    mc2 = MutCondition('61', 'Q', 'L')
    nras1 = Agent('NRAS', mutations=[mc1], db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', mutations=[mc2], db_refs={'HGNC': '7989'})
    nras3 = Agent('NRAS', mutations=[mc1], db_refs={'HGNC': 'dummy'})
    st1 = ActiveForm(nras1, 'gtpbound', True,
                     evidence=Evidence(text='foo'))
    st2 = ActiveForm(nras2, 'phosphatase', False,
                     evidence=Evidence(text='bar'))
    st3 = ActiveForm(nras3, 'gtpbound', False,
                     evidence=Evidence(text='baz'))
    assert st1.entities_match(st2)
    assert not st1.entities_match(st3)
    assert unicode_strs((st1, st2, st3))


def test_entities_match_gef():
    """Test matching of entities only, entities match on name and grounding."""
    sos1 = Agent('SOS1', db_refs={'HGNC': '123'})
    sos2 = Agent('SOS1', db_refs={'HGNC': '234'})
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': '7989'})
    st1 = Gef(sos1, nras1,
              evidence=Evidence(text='foo'))
    st2 = Gef(sos2, nras2,
              evidence=Evidence(text='bar'))
    st3 = Gef(sos1, nras2,
              evidence=Evidence(text='bar'))
    assert not st1.entities_match(st2), (st1.matches_key(), st2.matches_key())
    assert not st2.entities_match(st3)
    assert st1.entities_match(st3)
    assert unicode_strs((st1, st2, st3))


def test_entities_match_gap():
    """Test matching of entities only, entities match on name and grounding."""
    rasa1 = Agent('RASA1', db_refs={'HGNC': '123'})
    rasa2 = Agent('RASA1', db_refs={'HGNC': '234'})
    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs={'HGNC': 'dummy'})
    st1 = Gap(rasa1, nras1,
              evidence=Evidence(text='foo'))
    st2 = Gap(rasa2, nras2,
              evidence=Evidence(text='bar'))
    assert not st1.entities_match(st2)


def test_entities_match_complex():
    """Test matching of entities only, entities match on name and grounding."""
    ksr1 = Agent('KSR1', db_refs={'HGNC': '123'})
    ksr2 = Agent('KSR1', db_refs={'HGNC': '234'})
    braf1 = Agent('BRAF', db_refs={'HGNC': '345'})
    braf2 = Agent('BRAF', db_refs={'HGNC': '456'})
    map2k1 = Agent('MAP2K1', db_refs={'HGNC': '567'})
    map2k2 = Agent('MAP2K1', db_refs={'HGNC': '678'})
    st1 = Complex([ksr1, braf1, map2k1], evidence=Evidence(text='foo'))
    st2 = Complex([ksr2, braf2, map2k2], evidence=Evidence(text='bar'))
    st3 = Complex([braf2, map2k2, ksr2], evidence=Evidence(text='baz'))
    assert not st1.entities_match(st2)
    assert st2.entities_match(st3)
    assert not st3.entities_match(st1)


def test_agent_superfamily_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    nras = Agent('NRAS', db_refs={'HGNC': '7989'})
    assert nras.refinement_of(ras, hierarchies)
    assert not ras.refinement_of(nras, hierarchies)
    # The top-level list should contain only one statement, the gene-level
    # one, supported by the family one.


def test_agent_boundcondition_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    bc1 = BoundCondition(Agent('BRAF', db_refs={'HGNC': '1097'}), True)
    bc2 = BoundCondition(Agent('RAF1', db_refs={'HGNC': '9829'}), True)
    bc3 = BoundCondition(Agent('RAF1', db_refs={'HGNC': '9829'}), False)
    bc4 = BoundCondition(Agent('RAF', db_refs={'FPLX': 'RAF'}), True)

    nras1 = Agent('NRAS', db_refs={'HGNC': '7989'}, bound_conditions=[bc1])
    nras2 = Agent('NRAS', db_refs={'HGNC': '7989'}, bound_conditions=[bc2])
    nras3 = Agent('NRAS', db_refs={'HGNC': '7989'}, bound_conditions=[bc3])
    nras4 = Agent('NRAS', db_refs={'HGNC': '7989'})
    nras5 = Agent('NRAS', db_refs={'HGNC': '7989'},
                  bound_conditions=[bc4])

    # nras1 (bound to BRAF)
    assert not nras2.refinement_of(nras1, hierarchies)
    assert not nras3.refinement_of(nras1, hierarchies)
    assert not nras4.refinement_of(nras1, hierarchies)
    assert not nras5.refinement_of(nras1, hierarchies)
    # nras2 (bound to CRAF)
    assert not nras1.refinement_of(nras2, hierarchies)
    assert not nras3.refinement_of(nras2, hierarchies)  # Not bound condition
    assert not nras4.refinement_of(nras2, hierarchies)
    assert not nras5.refinement_of(nras2, hierarchies)
    # nras3 (not bound to CRAF)
    assert not nras1.refinement_of(nras3, hierarchies)
    assert not nras2.refinement_of(nras3, hierarchies)  # Not bound condition
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
    mek1 = Agent('MAP2K1', db_refs={'HGNC': 'asdf'},
                 mods=ModCondition('phosphorylation'))
    mek2 = Agent('MAP2K1', db_refs={'HGNC': 'asdf'},
                 mods=ModCondition('phosphorylation', position='218'))
    mek3 = Agent('MAP2K1', db_refs={'HGNC': 'asdf'},
                 mods=ModCondition('phosphorylation', position='222'))
    mek4 = Agent('MAP2K1', db_refs={'HGNC': 'asdf'},
                 mods=[ModCondition('phosphorylation', position='218'),
                       ModCondition('phosphorylation', position='222')])
    mek5 = Agent('MAP2K1', db_refs={'HGNC': 'asdf'},
                 mods=ModCondition('phosphorylation', 'serine', None))
    mek6 = Agent('MAP2K1', db_refs={'HGNC': 'asdf'},
                 mods=ModCondition('phosphorylation', 'serine', '218'))
    mek7 = Agent('MAP2K1', db_refs={'HGNC': 'asdf'},
                 mods=ModCondition('phosphorylation', 'serine', '222'))
    mek8 = Agent('MAP2K1', db_refs={'HGNC': 'asdf'},
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
    assert not mek3.refinement_of(mek2, hierarchies)  # Different site
    assert mek4.refinement_of(mek2, hierarchies)
    assert not mek5.refinement_of(mek2, hierarchies)  # Cross-relationship
    assert mek6.refinement_of(mek2, hierarchies)
    assert not mek7.refinement_of(mek2, hierarchies)  # Different site
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
    braf = Agent('BRAF', db_refs={'HGNC': 'braf'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': 'map2k1'})
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
    braf = Agent('BRAF', db_refs={'HGNC': 'braf'})
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
    raf = Agent('RAF', db_refs={'FPLX': 'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    mek = Agent('MEK', db_refs={'FPLX': 'MEK'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})

    st1 = Activation(raf, mek, 'kinase')
    st2 = Activation(braf, mek, 'kinase')
    st3 = Activation(raf, mek1, 'kinase')
    st4 = Activation(braf, mek1, 'kinase')
    st5 = Inhibition(braf, mek1, 'kinase')
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
    raf_k = Agent('RAF', activity=ActivityCondition('kinase', True),
                  db_refs={'FPLX': 'RAF'})
    raf_c = Agent('RAF', activity=ActivityCondition('catalytic', True),
                  db_refs={'FPLX': 'RAF'})
    raf_a = Agent('RAF', activity=ActivityCondition('activity', True),
                  db_refs={'FPLX': 'RAF'})
    mek = Agent('MEK', db_refs={'FPLX': 'MEK'})

    st1 = Activation(raf_k, mek, 'kinase')
    st2 = Inhibition(raf_k, mek, 'kinase')
    st3 = Activation(raf_c, mek, 'kinase')
    st4 = Activation(raf_k, mek, 'catalytic')
    st5 = Activation(raf_c, mek, 'activity')
    st6 = Activation(raf_a, mek, 'activity')

    assert not st1.refinement_of(st2, hierarchies)
    assert not st2.refinement_of(st1, hierarchies)
    assert st1.refinement_of(st3, hierarchies)
    assert st1.refinement_of(st4, hierarchies)
    assert st5.refinement_of(st6, hierarchies)
    assert st1.refinement_of(st6, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)
    assert not st4.refinement_of(st3, hierarchies)


def test_activitymod_refinement():
    mc1 = ModCondition('phosphorylation')
    mc2 = ModCondition('phosphorylation', 'S')
    mc3 = ModCondition('phosphorylation', 'S', '218')
    mc4 = ModCondition('phosphorylation', 'S', '222')
    mek_fam = Agent('MEK')
    mek1 = Agent('MAP2K1')
    p1 = ActiveForm(Agent('MEK', mods=[mc1], db_refs={'FPLX': 'MEK'}),
                    'kinase', True)
    p2 = ActiveForm(Agent('MEK', mods=[mc3], db_refs={'FPLX': 'MEK'}),
                    'kinase', True)
    p3 = ActiveForm(Agent('MAP2K1', mods=[mc1], db_refs={'HGNC': '6840'}),
                    'kinase', True)
    p4 = ActiveForm(Agent('MAP2K1', mods=[mc2], db_refs={'HGNC': '6840'}),
                    'kinase', True)
    p5 = ActiveForm(Agent('MAP2K1', mods=[mc3], db_refs={'HGNC': '6840'}),
                    'kinase', True)
    p6 = ActiveForm(Agent('MAP2K1', mods=[mc4], db_refs={'HGNC': '6840'}),
                    'kinase', True)
    p7 = ActiveForm(Agent('MAP2K1', mods=[mc3, mc4], db_refs={'HGNC': '6840'}),
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

    assert not p1.refinement_of(p2, hierarchies)
    assert p1.refinement_of(p3, hierarchies)
    assert p1.refinement_of(p4, hierarchies)
    assert p3.refinement_of(p4, hierarchies)
    assert not p4.refinement_of(p3, hierarchies)


def test_activatingsub_family_refinement():
    mc = MutCondition('12', 'G', 'D')
    ras = Agent('RAS', mutations=[mc], db_refs={'FPLX': 'RAS'})
    kras = Agent('KRAS', mutations=[mc], db_refs={'HGNC': '6407'})
    nras = Agent('NRAS', mutations=[mc], db_refs={'HGNC': '7989'})
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


def test_gef_family_refinement():
    sos = Agent('SOS', db_refs={'FPLX': 'SOS'})
    sos1 = Agent('SOS1', db_refs={'HGNC': '11187'})
    sos1_a = Agent('SOS1', activity=ActivityCondition('activity', True),
                   db_refs={'HGNC': '11187'})
    sos1_c = Agent('SOS1', activity=ActivityCondition('catalytic', True),
                   db_refs={'HGNC': '11187'})
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    kras = Agent('KRAS', db_refs={'HGNC': '6407'})
    # Statements
    st1 = Gef(sos, ras)
    st2 = Gef(sos1, ras)
    st3 = Gef(sos, kras)
    st4 = Gef(sos1, kras)
    st5 = Gef(sos1_a, kras)
    st6 = Gef(sos1_c, kras)
    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert st4.refinement_of(st1, hierarchies)
    assert st5.refinement_of(st1, hierarchies)
    assert st6.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert st4.refinement_of(st2, hierarchies)
    assert st5.refinement_of(st2, hierarchies)
    assert st6.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st4.refinement_of(st3, hierarchies)
    assert st5.refinement_of(st3, hierarchies)
    assert st6.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)
    assert st5.refinement_of(st4, hierarchies)
    assert st6.refinement_of(st4, hierarchies)
    # st5
    assert not st1.refinement_of(st5, hierarchies)
    assert not st2.refinement_of(st5, hierarchies)
    assert not st3.refinement_of(st5, hierarchies)
    assert not st4.refinement_of(st5, hierarchies)
    assert st6.refinement_of(st5, hierarchies)
    # st6
    assert not st5.refinement_of(st6, hierarchies)


def test_gap_family_refinement():
    rasa = Agent('RASA', db_refs={'FPLX': 'RASA'})
    rasa1 = Agent('RASA1', db_refs={'HGNC': '9871'})
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    kras = Agent('KRAS', db_refs={'HGNC': '6407'})
    rasa1_a = Agent('RASA1', activity=ActivityCondition('activity', True),
                    db_refs={'HGNC': '9871'})
    rasa1_c = Agent('RASA1', activity=ActivityCondition('catalytic', True),
                    db_refs={'HGNC': '9871'})
    # Statements
    st1 = Gap(rasa, ras)
    st2 = Gap(rasa1, ras)
    st3 = Gap(rasa, kras)
    st4 = Gap(rasa1, kras)
    st5 = Gap(rasa1_a, kras)
    st6 = Gap(rasa1_c, kras)
    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert st4.refinement_of(st1, hierarchies)
    assert st5.refinement_of(st1, hierarchies)
    assert st6.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert st4.refinement_of(st2, hierarchies)
    assert st5.refinement_of(st2, hierarchies)
    assert st6.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st4.refinement_of(st3, hierarchies)
    assert st5.refinement_of(st3, hierarchies)
    assert st6.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)
    assert st5.refinement_of(st4, hierarchies)
    assert st6.refinement_of(st4, hierarchies)
    # st5
    assert not st1.refinement_of(st5, hierarchies)
    assert not st2.refinement_of(st5, hierarchies)
    assert not st3.refinement_of(st5, hierarchies)
    assert not st4.refinement_of(st5, hierarchies)
    assert st6.refinement_of(st5, hierarchies)
    # st6
    assert not st5.refinement_of(st6, hierarchies)


def test_complex_family_refinement():
    raf = Agent('RAF', db_refs={'FPLX': 'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    raf1 = Agent('RAF1', db_refs={'HGNC': '9829'})
    mek = Agent('MEK', db_refs={'FPLX': 'MEK'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})

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


def test_related_complex_refinement():
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    kras = Agent('KRAS', db_refs={'HGNC': '6407'})
    hras = Agent('HRAS', db_refs={'HGNC': '5173'})
    st1 = Complex([kras, hras])
    st2 = Complex([kras, ras])
    st3 = Complex([hras, kras])
    st4 = Complex([ras, kras])
    assert st1.refinement_of(st2, hierarchies)
    assert st3.refinement_of(st4, hierarchies)
    assert st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st1, hierarchies)
    assert not st4.refinement_of(st1, hierarchies)


def test_big_complex_refinement():
    c1js = ('{"type": "Complex", "members": [{"name": "MMP2", "db_refs": '
            '{"UP": "P08253", "TEXT": "MMP-2", "HGNC": "7166"}},'
            '{"name": "MMP9", "db_refs": {"UP": "P14780", "TEXT": "MMP-9", '
            '"HGNC": "7176"}}, {"name": "ERVK-6", "db_refs": {"UP": "Q9Y6I0", '
            '"HGNC": "13915"}}, {"name": "ERVK-7", "db_refs": {"UP": "P63131", '
            '"HGNC": "31828"}}, {"name": "ERVK-24", "db_refs": '
            '{"UP": "P63129", "HGNC": "39038"}}, {"name": "ERVK-9", "db_refs": '
            '{"UP": "P63127", "HGNC": "39005"}}, {"name": "ERVK-25", '
            '"db_refs": {"UP": "P63125", "HGNC": "39039"}}, '
            '{"name": "HERV-K104", "db_refs": {"UP": "P63124"}}, '
            '{"name": "ERVK-18", "db_refs": {"UP": "P63123", '
            '"HGNC": "39025"}}, {"name": "ERVK-8", "db_refs": {"UP": "P63122", '
            '"HGNC": "32302"}}, {"name": "HERVK_113", "db_refs": '
            '{"UP": "P63121"}}, {"name": "ERVK-19", "db_refs": '
            '{"UP": "P63120", "HGNC": "39026"}}, {"name": "F2RL1", "db_refs": '
            '{"UP": "P55085", "HGNC": "3538"}}, {"name": "ERVK-10", "db_refs": '
            '{"UP": "P10265", "HGNC": "39004"}}], "belief": 1.0, "evidence": '
            '[{"source_api": "sparser", "text": "The results suggest that for '
            'the completion of invasion and migration of H-ras MCF10A cells, '
            'not only matrix-degrading proteinase activities of MMP-2 and '
            'MMP-9 but also other intracellular events through activation of '
            'ERKs signaling pathways are essential.", "annotations": '
            '{"found_by": "BIO-ACTIVITY"}}], "id": '
            '"c41869d1-d48b-4756-8a76-052515ae3a6d"}')
    c2js = ('{"type": "Complex", "members": [{"name": "MMP3", "db_refs": '
            '{"UP": "P08254", "TEXT": "MMP-3", "HGNC": "7173"}}, '
            '{"name": "MMP7", "db_refs": {"UP": "P09237", "TEXT": "MMP-7", '
            '"HGNC": "7174"}}, {"name": "ERVK-6", "db_refs": {"UP": "Q9Y6I0", '
            '"HGNC": "13915"}}, {"name": "ERVK-7", "db_refs": {"UP": "P63131", '
            '"HGNC": "31828"}}, {"name": "ERVK-24", "db_refs": '
            '{"UP": "P63129", "HGNC": "39038"}}, {"name": "ERVK-9", "db_refs": '
            '{"UP": "P63127", "HGNC": "39005"}}, {"name": "ERVK-25", '
            '"db_refs": {"UP": "P63125", "HGNC": "39039"}}, '
            '{"name": "HERV-K104", "db_refs": {"UP": "P63124"}}, '
            '{"name": "ERVK-18", "db_refs": {"UP": "P63123", '
            '"HGNC": "39025"}}, {"name": "ERVK-8", "db_refs": {"UP": "P63122", '
            '"HGNC": "32302"}}, {"name": "HERVK_113", "db_refs": '
            '{"UP": "P63121"}}, {"name": "ERVK-19", "db_refs": {"UP": "P63120",'
            ' "HGNC": "39026"}}, {"name": "F2RL1", "db_refs": {"UP": "P55085", '
            '"HGNC": "3538"}}, {"name": "ERVK-10", "db_refs": {"UP": "P10265", '
            '"HGNC": "39004"}}], "belief": 1.0, "evidence": [{"source_api": '
            '"sparser", "text": "Proteinase activities of MMP-3 and MMP-7 were '
            'quantified by means of casein zymography (,  and ).", '
            '"annotations": {"found_by": "BIO-ACTIVITY"}}], "id": '
            '"5a03d72c-755e-4262-b8c8-9f632305c703"}')
    c1 = Statement._from_json(json.loads(c1js))
    c2 = Statement._from_json(json.loads(c2js))
    c1.refinement_of(c2, hierarchies)


def test_mtor_rictor_refinement():
    mt = Agent('MTOR')
    mtm = Agent('MTOR', mods=[ModCondition('phosphorylation')])
    rt = Agent('RICTOR')
    c1 = Complex([mt, rt, rt])
    c2 = Complex([mtm, mt, mt])
    assert not c1.refinement_of(c2, hierarchies)
    assert not c2.refinement_of(c1, hierarchies)


@raises(InvalidResidueError)
def test_residue_mod_condition():
    ModCondition('phosphorylation', 'xyz')


@raises(InvalidResidueError)
def test_residue_mod():
    Phosphorylation(Agent('a'), Agent('b'), 'xyz')


@raises(InvalidResidueError)
def test_residue_selfmod():
    Autophosphorylation(Agent('a'), 'xyz')


def test_valid_mod_residue():
    mc = ModCondition('phosphorylation', 'serine')
    assert mc.residue == 'S'
    assert unicode_strs(mc)


def test_valid_residue():
    assert get_valid_residue('serine') == 'S'
    assert get_valid_residue('ser') == 'S'
    assert get_valid_residue('Serine') == 'S'
    assert get_valid_residue('SERINE') == 'S'


def test_modcondition_order_actmod():
    mc1 = ModCondition('phoshporylation', 'S', '222')
    mc2 = ModCondition('phoshporylation', 'S', '224')
    p1 = ActiveForm(Agent('MAP2K1', mods=[mc1, mc2]), 'kinase', True)
    p2 = ActiveForm(Agent('MAP2K1', mods=[mc2, mc1]), 'kinase', True)
    assert p1.matches(p2)
    assert unicode_strs((p1, p2))


def test_modcondition_order_agent():
    mc1 = ModCondition('phoshporylation', 'S', '222')
    mc2 = ModCondition('phoshporylation', 'S', '224')
    p1 = Agent('MAP2K1', mods=[mc1, mc2])
    p2 = Agent('MAP2K1', mods=[mc2, mc1])
    assert p1.matches(p2)
    assert unicode_strs((p1, p2))


def test_eq_mut():
    assert MutCondition('600', 'V', 'E').equals(MutCondition('600', 'V', 'E'))
    assert not MutCondition('600', 'V', 'E').equals(
        MutCondition('600', 'V', 'D'))
    return


def test_mut_hgvs():
    mc = MutCondition('600', 'V', 'E')
    assert mc.to_hgvs() == 'p.Val600Glu'


def test_eq_agent():
    assert Agent('one').equals(Agent('one'))
    assert not Agent('one').equals(Agent('two'))
    assert not Agent('one', db_refs={'UP': '123'}).equals(
           Agent('one', db_refs={'UP': '999'}))
    assert Agent('one', mods=[ModCondition('phosphorylation')]).equals(
           Agent('one', mods=[ModCondition('phosphorylation')]))
    assert not Agent('one', mods=[ModCondition('phosphorylation')]).equals(
           Agent('one', mods=[ModCondition('ubiquitination')]))
    assert Agent('one', mutations=[MutCondition('600', 'V', 'E')]).equals(
           Agent('one', mutations=[MutCondition('600', 'V', 'E')]))
    assert not Agent('one', mutations=[MutCondition('600', 'V', 'E')]).equals(
           Agent('one', mutations=[MutCondition('600', 'V', 'D')]))
    assert Agent('one',
                 bound_conditions=[BoundCondition(Agent('two'), True)]).equals(
           Agent('one',
                 bound_conditions=[BoundCondition(Agent('two'), True)]))
    assert not Agent('one',
                     bound_conditions=[BoundCondition(Agent('two'),
                                                      True)]).equals(
           Agent('one', bound_conditions=[BoundCondition(Agent('two'),
                                                         False)]))
    assert not Agent('one', bound_conditions=[BoundCondition(Agent('two'),
                                                             True)]).equals(
           Agent('one', bound_conditions=[BoundCondition(Agent('three'),
                                                         True)]))
    return


def test_eq_stmt():
    ev1 = Evidence(text='1')
    ev2 = Evidence(text='2')
    assert Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
            Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]))
    assert not Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
            Phosphorylation(Agent('a'), Agent('b'), evidence=[ev2]))
    assert not Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
            Phosphorylation(Agent('a'), Agent('c'), evidence=[ev2]))
    assert not Phosphorylation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
            Phosphorylation(Agent('a'), Agent('b'), 'S', evidence=[ev2]))
    assert Complex([Agent('a'), Agent('b')], evidence=[ev1]).equals(
           Complex([Agent('a'), Agent('b')], evidence=[ev1]))
    assert not Complex([Agent('a'), Agent('b')], evidence=[ev1]).equals(
           Complex([Agent('a'), Agent('b')], evidence=[ev2]))
    assert Activation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
           Activation(Agent('a'), Agent('b'), evidence=[ev1]))
    assert not Activation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
           Activation(Agent('a'), Agent('c'), evidence=[ev1]))
    assert not Activation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
           Activation(Agent('a'), Agent('b'), 'kinase', evidence=[ev1]))
    assert not Activation(Agent('a'), Agent('b'), evidence=[ev1]).equals(
           Activation(Agent('a'), Agent('b'), evidence=[ev2]))


def test_serialize():
    ev1 = Evidence(text='1\U0001F4A9')
    st = Phosphorylation(Agent('a\U0001F4A9'), Agent('b'), evidence=[ev1])
    jstr = st.to_json()
    st2 = stmts_from_json([jstr])[0]
    assert st.equals(st2)
    assert unicode_strs((ev1, st, st2))
    assert st.evidence[0].source_hash == st2.evidence[0].source_hash


def test_location_refinement():
    a1 = Agent('a', location='plasma membrane')
    a2 = Agent('a', location='cell')
    a3 = Agent('a', location='cytoplasm')
    a4 = Agent('a')
    a5 = Agent('a')

    assert a1.refinement_of(a2, hierarchies)
    assert not a2.refinement_of(a3, hierarchies)
    assert a4.refinement_of(a5, hierarchies)
    assert not a1.refinement_of(a3, hierarchies)
    assert not a3.refinement_of(a1, hierarchies)
    assert a2.refinement_of(a4, hierarchies)
    assert a3.refinement_of(a4, hierarchies)


def test_activity_refinement():
    a1 = Agent('a', activity=ActivityCondition('kinase', True))
    a2 = Agent('a', activity=ActivityCondition('activity', True))
    a3 = Agent('a', activity=ActivityCondition('catalytic', True))
    a4 = Agent('a')
    a5 = Agent('a', activity=ActivityCondition('catalytic', False))
    a6 = Agent('a', activity=ActivityCondition('kinase', False))

    assert a1.refinement_of(a2, hierarchies)
    assert not a2.refinement_of(a3, hierarchies)
    assert not a4.refinement_of(a1, hierarchies)
    assert a1.refinement_of(a3, hierarchies)
    assert a3.refinement_of(a2, hierarchies)
    assert not a3.refinement_of(a1, hierarchies)
    assert a1.refinement_of(a4, hierarchies)
    assert a2.refinement_of(a4, hierarchies)
    assert a5.refinement_of(a4, hierarchies)
    assert not a5.refinement_of(a3, hierarchies)
    assert not a5.refinement_of(a1, hierarchies)
    assert a6.refinement_of(a5, hierarchies)
    assert not a5.refinement_of(a6, hierarchies)


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
    assert st3.refinement_of(st8, hierarchies)
    assert st1.refinement_of(st2, hierarchies)
    assert st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert st4.refinement_of(st5, hierarchies)
    assert st6.refinement_of(st5, hierarchies)
    assert not st1.refinement_of(st7, hierarchies)
    assert st7.refinement_of(st4, hierarchies)
    assert st8.refinement_of(st9, hierarchies)
    assert st7.refinement_of(st9, hierarchies)


def test_decrease_amt_refinement():
    raf = Agent('RAF', db_refs={'FPLX': 'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    brafk = Agent('BRAF', activity=ActivityCondition('kinase', True),
                  db_refs={'HGNC': '1097'})
    raf1 = Agent('RAF1', db_refs={'HGNC': '9829'})
    mek = Agent('MEK', db_refs={'FPLX': 'MEK'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})

    st1 = DecreaseAmount(raf, mek)
    st2 = DecreaseAmount(braf, mek)
    st3 = DecreaseAmount(raf, mek1)
    st4 = DecreaseAmount(brafk, mek1)

    assert unicode_strs((st1, st2, st3, st4))
    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert st4.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert st4.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st4.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)


def test_increase_amt_refinement():
    raf = Agent('RAF', db_refs={'FPLX': 'RAF'})
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    brafk = Agent('BRAF', activity=ActivityCondition('kinase', True),
                  db_refs={'HGNC': '1097'})
    raf1 = Agent('RAF1', db_refs={'HGNC': '9829'})
    mek = Agent('MEK', db_refs={'FPLX': 'MEK'})
    mek1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})

    st1 = IncreaseAmount(raf, mek)
    st2 = IncreaseAmount(braf, mek)
    st3 = IncreaseAmount(raf, mek1)
    st4 = IncreaseAmount(brafk, mek1)

    assert unicode_strs((st1, st2, st3, st4))
    # st1
    assert st2.refinement_of(st1, hierarchies)
    assert st3.refinement_of(st1, hierarchies)
    assert st4.refinement_of(st1, hierarchies)
    # st2
    assert not st1.refinement_of(st2, hierarchies)
    assert not st3.refinement_of(st2, hierarchies)
    assert st4.refinement_of(st2, hierarchies)
    # st3
    assert not st1.refinement_of(st3, hierarchies)
    assert not st2.refinement_of(st3, hierarchies)
    assert st4.refinement_of(st3, hierarchies)
    # st4
    assert not st1.refinement_of(st4, hierarchies)
    assert not st2.refinement_of(st4, hierarchies)
    assert not st3.refinement_of(st4, hierarchies)


def test_complex_refinement_order():
    st1 = Complex([Agent('MED23'), Agent('ELK1')])
    st2 = Complex([Agent('ELK1', mods=[ModCondition('phosphorylation')]),
                   Agent('MED23')])
    assert st2.refinement_of(st1, hierarchies)
    assert not st1.refinement_of(st2, hierarchies)


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


def test_mod_condition_is_mod():
    mc1 = ModCondition('ubiquitination', 'K', '99', True)
    mc2 = ModCondition('ubiquitination', 'K', '99', False)
    assert not mc1.refinement_of(mc2, hierarchies)


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

    st = Activation(ag, ag, 'activity', evidence=ev)
    print(st)
    print(repr(st))

    st = Inhibition(ag, ag, 'activity', evidence=ev)
    print(st)
    print(repr(st))

    st = ActiveForm(ag, 'activity', True)
    print(st)
    print(repr(st))

    st = HasActivity(ag, 'activity', True)
    print(st)
    print(repr(st))

    st = Gef(ag, ag, evidence=ev)
    print(st)
    print(repr(st))

    st = Gap(ag, ag, evidence=ev)
    print(st)
    print(repr(st))

    st = Complex([ag, ag], evidence=ev)
    print(st)
    print(repr(st))


def test_modtype_to_modclass():
    cls = modtype_to_modclass.get('farnesylation')
    assert cls == Farnesylation
    cls = modtype_to_modclass.get('deubiquitination')
    assert cls == Deubiquitination


def test_modclass_to_modtype():
    modtype = modclass_to_modtype.get(Depalmitoylation)
    assert modtype == 'depalmitoylation'
    modtype = modclass_to_modtype.get(Phosphorylation)
    assert modtype == 'phosphorylation'


def test_modtype_to_inverse():
    modtype_inv = modtype_to_inverse.get('ubiquitination')
    assert modtype_inv == 'deubiquitination'
    modtype_inv = modtype_to_inverse.get('dephosphorylation')
    assert modtype_inv == 'phosphorylation'


def test_mut_refinement():
    mc1 = MutCondition('600', 'V', 'E')
    mc2 = MutCondition('600', 'V', None)
    mc3 = MutCondition(None, 'V', 'E')
    mc4 = MutCondition(None, None, None)
    assert mc1.refinement_of(mc2)
    assert mc1.refinement_of(mc3)
    assert mc1.refinement_of(mc3)
    assert not mc2.refinement_of(mc1)
    assert not mc2.refinement_of(mc3)
    assert mc2.refinement_of(mc4)
    assert not mc3.refinement_of(mc2)
    assert mc3.refinement_of(mc4)
    assert not mc4.refinement_of(mc1)
    assert not mc4.refinement_of(mc2)
    assert not mc4.refinement_of(mc3)


def test_mut_agent_refinement():
    mc1 = MutCondition('600', 'V', 'E')
    mc2 = MutCondition('600', 'V', None)
    mc3 = MutCondition(None, 'V', 'E')
    mc4 = MutCondition(None, None, None)
    a1 = Agent('a', mutations=[mc1])
    a2 = Agent('a', mutations=[mc4])
    assert a1.refinement_of(a2, hierarchies)
    assert not a2.refinement_of(a1, hierarchies)


def test_conversion_init():
    Conversion(Agent('RAS'), Agent('GTP'), Agent('GDP'))


def test_conversion_refinement():
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    hras = Agent('HRAS', db_refs={'HGNC': '5173'})
    gtp = Agent('GTP')
    gdp = Agent('GDP')
    st1 = Conversion(ras, gtp, gdp)
    st2 = Conversion(hras, gtp, gdp)
    assert st2.refinement_of(st1, hierarchies)
    assert not st1.refinement_of(st2, hierarchies)


def test_conversion_set_agent_list():
    ag = [Agent('A%d' % i) for i in range(5)]
    st1 = Conversion(None, [ag[0]], [ag[1]])
    st2 = Conversion(ag[0], [ag[1], ag[2]], [ag[3]])
    st3 = Conversion(ag[0], [ag[1]], [ag[2], ag[3]])
    st1.set_agent_list([None] + ag[:2])
    st2.set_agent_list(ag[:4])
    st3.set_agent_list(ag[:4])
    assert st1.subj is None
    assert len(st1.obj_from) == 1
    assert len(st1.obj_to) == 1
    assert len(st2.obj_from) == 2
    assert len(st2.obj_to) == 1
    assert len(st3.obj_from) == 1
    assert len(st3.obj_to) == 2


def test_get_act_condition():
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Activation(braf, mek)
    ac = stmt._get_activity_condition()
    assert ac.activity_type == 'activity'
    assert ac.is_active
    stmt = Inhibition(braf, mek, 'kinase')
    ac = stmt._get_activity_condition()
    assert ac.activity_type == 'kinase'
    assert not ac.is_active


def test_bound_condition_matches():
    bcs = [BoundCondition(Agent('a'), True),
           BoundCondition(Agent('b'), True),
           BoundCondition(Agent('a'), False),
           BoundCondition(Agent('a'), True)]
    assert not bcs[0].matches(bcs[1])
    assert not bcs[0].matches(bcs[2])
    assert bcs[0].matches(bcs[3])


def test_influence_polarity():
    st = Influence(Event(Concept('a')), Event(Concept('b')))
    assert st.overall_polarity() is None
    st.subj.delta = QualitativeDelta(polarity=None, adjectives=None)
    assert st.overall_polarity() is None
    st.obj.delta = QualitativeDelta(polarity=None, adjectives=None)
    assert st.overall_polarity() is None
    st.subj.delta.set_polarity(1)
    assert st.overall_polarity() == 1
    st.subj.delta.set_polarity(-1)
    assert st.overall_polarity() == -1
    st.obj.delta.set_polarity(1)
    assert st.overall_polarity() == -1
    st.obj.delta.set_polarity(-1)
    assert st.overall_polarity() == 1
    st.subj.delta.set_polarity(1)
    assert st.overall_polarity() == -1
    st.obj.delta.set_polarity(1)
    assert st.overall_polarity() == 1, st


def test_association_polarity():
    st = Association([Event(Concept('a')), Event(Concept('b'))])
    assert st.overall_polarity() is None
    st.members[0].delta.set_polarity(1)
    assert st.overall_polarity() == 1
    st.members[1].delta.set_polarity(1)
    assert st.overall_polarity() == 1
    st.members[0].delta.set_polarity(-1)
    assert st.overall_polarity() == -1
    st.members[1].delta.set_polarity(-1)
    assert st.overall_polarity() == 1
    st.members[0].delta = QualitativeDelta(polarity=None, adjectives=None)
    assert st.overall_polarity() == -1
    st.members[1].delta = QualitativeDelta(polarity=None, adjectives=None)
    assert st.overall_polarity() is None


def test_concept_init():
    c1 = Concept('x')
    assert c1.name == 'x'
    c2 = Concept('y', db_refs={'TEXT': 'y'})
    assert c2.name == 'y'
    assert c2.db_refs['TEXT'] == 'y'


def test_concept_matches():
    # matches
    assert Concept('x').matches(Concept('x', db_refs={'TEXT': 'x'}))
    assert not Concept('x').matches(Concept('y'))
    assert Concept('x', db_refs={'EIDOS': 'x'}).matches(
        Concept('y', db_refs={'EIDOS': 'x'}))
    # entity_matches
    assert Concept('x').entity_matches(Concept('x', db_refs={'TEXT': 'x'}))
    assert not Concept('x').entity_matches(Concept('y'))
    # equals
    assert Concept('x').equals(Concept('x'))
    assert not Concept('x').equals(Concept('y'))
    assert Concept('x', db_refs={'TEXT': 'x'}).equals(
           Concept('x', db_refs={'TEXT': 'x'}))
    assert not Concept('x').equals(Concept('x', db_refs={'TEXT': 'x'}))


def test_concept_get_grounding():
    d1 = {'TEXT': 'a'}
    d2 = {'TEXT': 'b', 'UN': 'c'}
    d3 = {'TEXT': 'x', 'UN': 'y', 'HUME': 'z'}
    d4 = {'TEXT': 'b', 'HUME': 'a'}
    d5 = {'UN': [('a', 1.0), ('b', 0.8)]}
    d6 = {'UN': [('b', 0.8), ('a', 1.0)]}
    d7 = {'UN': []}
    d8 = {'HUME': [('a', 1.0), ('b', 0.8)]}
    assert Concept('a', db_refs=d1).get_grounding() == (None, None)
    assert Concept('b', db_refs=d2).get_grounding() == ('UN', 'c')
    assert Concept('c', db_refs=d3).get_grounding() == ('UN', 'y')
    assert Concept('d', db_refs=d4).get_grounding() == ('HUME', 'a')
    assert Concept('e', db_refs=d5).get_grounding() == ('UN', 'a')
    assert Concept('f', db_refs=d6).get_grounding() == ('UN', 'a')
    assert Concept('g', db_refs=d7).get_grounding() == (None, None)
    assert Concept('h', db_refs=d8).get_grounding() == ('HUME', 'a')


def test_concept_isa_eid():
    eidos_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/eidos/eidos_ontology.rdf')
    hm = HierarchyManager(eidos_ont, True, True)
    c1 = Concept('a', db_refs={'UN': [('UN/events/human/conflict', 1.0)]})
    c2 = Concept('b', db_refs={'UN': [('UN/events/human', 1.0)]})
    print(c1.get_grounding())
    print(c2.get_grounding())
    assert c1.refinement_of(c2, {'entity': hm})
    assert not c2.refinement_of(c1, {'entity': hm})


def test_concept_opposite_eid():
    eidos_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/eidos/eidos_ontology.rdf')
    hm = HierarchyManager(eidos_ont, True, True)
    c1 = Concept('a', db_refs={'UN':
                               [('UN/entities/human/food/food_insecurity',
                                 1.0)]})
    c2 = Concept('b', db_refs={'UN':
                               [('UN/entities/human/food/food_security',
                                 1.0)]})
    assert c1.is_opposite(c2, {'entity': hm})
    assert c2.is_opposite(c1, {'entity': hm})


def test_concept_isa_cwms():
    trips_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/cwms/trips_ontology.rdf')
    hm = HierarchyManager(trips_ont, True, True)
    c1 = Concept('a', db_refs={'CWMS': 'ONT::TRUCK'})
    c2 = Concept('b', db_refs={'CWMS': 'ONT::VEHICLE'})
    assert c1.refinement_of(c2, {'entity': hm})
    assert not c2.refinement_of(c1, {'entity': hm})


def test_concept_isa_hume():
    hume_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            '../sources/hume/hume_ontology.rdf')
    hm = HierarchyManager(hume_ont, True, True)
    c1 = Concept('a',
                 db_refs={'HUME': 'entity/rule/law'})
    c2 = Concept('b', db_refs={'HUME': 'entity/rule'})
    assert c1.refinement_of(c2, {'entity': hm})
    assert not c2.refinement_of(c1, {'entity': hm})


def test_influence_refinement_of():
    influence_association_refinement_of('Influence')


def test_association_refinement_of():
    influence_association_refinement_of('Association')


def influence_association_refinement_of(stmt_type):
    c1 = Event(Concept('production'))
    c2 = Event(Concept('price'))
    nopol_noadj = QualitativeDelta(polarity=None, adjectives=None)
    nopol_adj = QualitativeDelta(polarity=None, adjectives=['small'])
    neg_noadj = QualitativeDelta(polarity=-1, adjectives=None)
    neg_adj = QualitativeDelta(polarity=-1, adjectives=['small'])
    pos_noadj = QualitativeDelta(polarity=1, adjectives=None)
    pos_adj = QualitativeDelta(polarity=1, adjectives=['small'])
    pos_adj2 = QualitativeDelta(polarity=1, adjectives=['small', 'tiny'])
    pos_adj3 = QualitativeDelta(polarity=1, adjectives=['significant', 'large'])

    def S(x, y, stmt_type):
        cc1 = deepcopy(c1)
        cc2 = deepcopy(c2)
        cc1.delta = x
        cc2.delta = y
        if stmt_type == 'Influence':
            return Influence(cc1, cc2)
        elif stmt_type == 'Association':
            return Association([cc1, cc2])

    # Has polarity vs doesn't have polarity
    assert S(pos_adj, nopol_noadj, stmt_type).refinement_of(
        S(nopol_adj, nopol_noadj, stmt_type), hierarchies)
    assert S(nopol_adj, pos_noadj, stmt_type).refinement_of(
        S(nopol_adj, nopol_noadj, stmt_type), hierarchies)
    assert S(nopol_adj, neg_noadj, stmt_type).refinement_of(
        S(nopol_adj, nopol_noadj, stmt_type), hierarchies)
    # Has adjective vs doesn't have adjective
    assert S(neg_adj, nopol_adj, stmt_type).refinement_of(
        S(neg_noadj, nopol_adj, stmt_type), hierarchies)
    assert S(pos_adj, nopol_adj, stmt_type).refinement_of(
        S(pos_noadj, nopol_adj, stmt_type), hierarchies)
    assert S(nopol_adj, nopol_adj, stmt_type).refinement_of(
        S(nopol_noadj, nopol_adj, stmt_type), hierarchies)
    # Opposite polarity
    assert not S(neg_noadj, nopol_noadj, stmt_type).refinement_of(
            S(pos_noadj, nopol_noadj, stmt_type), hierarchies)
    # Adjectives can be in any relation
    assert S(pos_adj2, nopol_noadj, stmt_type).refinement_of(
        S(pos_adj, nopol_noadj, stmt_type), hierarchies)
    assert S(pos_adj3, nopol_noadj, stmt_type).refinement_of(
        S(pos_adj, nopol_noadj, stmt_type), hierarchies)
    assert S(pos_adj2, nopol_noadj, stmt_type).refinement_of(
        S(pos_adj3, nopol_noadj, stmt_type), hierarchies)

    # Equivalences
    assert not S(nopol_adj, neg_noadj, stmt_type).equals(
            S(nopol_adj, nopol_noadj, stmt_type))
    assert not S(nopol_adj, neg_noadj, stmt_type).equals(
            S(nopol_noadj, neg_noadj, stmt_type))
    pos_adj4 = QualitativeDelta(polarity=1, adjectives=['large', 'significant'])
    assert S(pos_adj3, nopol_noadj, stmt_type).equals(
        S(pos_adj4, nopol_noadj, stmt_type))

    # Matches keys
    assert S(nopol_adj, neg_noadj, stmt_type).matches_key() != \
        S(nopol_adj, nopol_noadj, stmt_type).matches_key()
    assert S(nopol_adj, neg_noadj, stmt_type).matches_key() == \
        S(nopol_noadj, neg_noadj, stmt_type).matches_key()
    pos_adj4 = QualitativeDelta(polarity=1, adjectives=['large', 'significant'])
    assert S(pos_adj3, nopol_noadj, stmt_type).matches_key() == \
        S(pos_adj4, nopol_noadj, stmt_type).matches_key()
    assert S(pos_adj, neg_adj, stmt_type).matches_key() == \
        S(neg_adj, pos_adj, stmt_type).matches_key()
    assert S(pos_adj, pos_adj, stmt_type).matches_key() == \
        S(neg_adj, neg_adj, stmt_type).matches_key()

    # Contradicts
    assert S(pos_adj, neg_noadj, stmt_type).contradicts(
        S(neg_adj, neg_noadj, stmt_type), hierarchies)
    assert S(pos_adj, neg_noadj, stmt_type).contradicts(
        S(pos_adj, pos_noadj, stmt_type), hierarchies)
    assert not S(pos_adj, neg_noadj, stmt_type).contradicts(
            S(pos_adj, neg_adj, stmt_type), hierarchies)
    assert not S(pos_adj, neg_noadj, stmt_type).contradicts(
            S(neg_adj, pos_adj, stmt_type), hierarchies)


def test_association_contradicts():
    eidos_ont = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '../sources/eidos/eidos_ontology.rdf')
    hm = HierarchyManager(eidos_ont, True, True)
    hierarchies = {'entity': hm}
    sn = Event(Concept('food security',
                       db_refs={'UN': 'UN/entities/human/food/food_security'}),
               delta=QualitativeDelta(polarity=-1))
    sp = Event(Concept('food security',
                       db_refs={'UN': 'UN/entities/human/food/food_security'}),
               delta=QualitativeDelta(polarity=1))
    ip = Event(Concept('food insecurity',
                       db_refs={'UN':
                                'UN/entities/human/food/food_insecurity'}),
               delta=QualitativeDelta(polarity=1))
    prp = Event(Concept('production'), delta=QualitativeDelta(polarity=1))
    prn = Event(Concept('production'), delta=QualitativeDelta(polarity=-1))
    assert Association([sn, prp]).contradicts(Association([sn, prn]),
                                              hierarchies)
    assert Association([prp, sn]).contradicts(Association([sn, prn]),
                                              hierarchies)
    assert Association([prn, sp]).contradicts(Association([sn, prn]),
                                              hierarchies)
    assert Association([sn, prp]).contradicts(Association([ip, prn]),
                                              hierarchies)
    assert Association([sn, sp]).contradicts(Association([ip, sn]),
                                             hierarchies)
    assert Association([ip, sp]).contradicts(Association([sp, sp]),
                                             hierarchies)
    assert Association([ip, sp]).contradicts(Association([sn, sn]),
                                             hierarchies)


def test_modification_contradicts():
    st1 = Phosphorylation(Agent('a'), Agent('b'))
    st2 = Dephosphorylation(Agent('a'), Agent('b'))
    st3 = Dephosphorylation(Agent('a'), Agent('b'), 'S')
    st4 = Dephosphorylation(Agent('a'), Agent('b'), 'S', '123')
    st5 = Phosphorylation(Agent('a'), Agent('b'), 'S', '123')
    st6 = Phosphorylation(Agent('a'), Agent('b'), 'T', '234')
    st7 = Phosphorylation(None, Agent('b'))
    st8 = Dephosphorylation(Agent('a'), Agent('c'))

    assert st1.contradicts(st2, hierarchies)
    assert not st1.contradicts(st7, hierarchies)
    assert not st1.contradicts(st5, hierarchies)
    assert not st4.contradicts(st6, hierarchies)
    assert st4.contradicts(st5, hierarchies)
    assert not st3.contradicts(st6, hierarchies)
    assert not st1.contradicts(st8, hierarchies)
    # TODO: add tests with Agent refinement
    # TODO: add tests with other Modification types


def test_regulate_amount_contradicts():
    st1 = IncreaseAmount(Agent('a'), Agent('b'))
    st2 = DecreaseAmount(Agent('a'), Agent('b'))
    st3 = DecreaseAmount(Agent('a'), Agent('c'))
    st4 = IncreaseAmount(Agent('b'), Agent('a'))
    assert st1.contradicts(st2, hierarchies)
    assert not st1.contradicts(st3, hierarchies)
    assert not st1.contradicts(st4, hierarchies)
    assert not st2.contradicts(st4, hierarchies)


def test_regulate_activity_contradicts():
    st1 = Activation(Agent('a'), Agent('b'))
    st2 = Inhibition(Agent('a'), Agent('b'))
    st3 = Inhibition(Agent('a'), Agent('c'))
    st4 = Activation(Agent('b'), Agent('a'))
    assert st1.contradicts(st2, hierarchies)
    assert not st1.contradicts(st3, hierarchies)
    assert not st1.contradicts(st4, hierarchies)
    assert not st2.contradicts(st4, hierarchies)


def test_active_form_contradicts():
    mc1 = ModCondition('phosphorylation', 'S', None, True)
    mc2 = ModCondition('phosphorylation', 'S', None, False)
    st1 = ActiveForm(Agent('a', mods=[mc1]), 'kinase', True)
    st2 = ActiveForm(Agent('a', mods=[mc1]), 'kinase', False)
    st3 = ActiveForm(Agent('a', mods=[mc1]), 'activity', False)
    st4 = ActiveForm(Agent('a', mods=[mc2]), 'activity', False)
    assert st1.contradicts(st2, hierarchies)
    assert not st1.contradicts(st3, hierarchies)
    assert not st3.contradicts(st4, hierarchies)
    assert not st1.contradicts(st4, hierarchies)


def test_agent_list_with_bound_condition_agents():
    eg = Agent('EGFR', bound_conditions=[BoundCondition(Agent('EGF'), True)])
    stmt = Phosphorylation(None, eg)
    agents = stmt.agent_list_with_bound_condition_agents()
    assert agents[0] is None
    assert agents[1].name == 'EGFR'
    assert agents[2].name == 'EGF'


def test_context_bool_equal():
    assert not BioContext()
    assert BioContext(cell_type=RefContext(name='x'))
    assert not WorldContext()
    assert WorldContext(time=TimeContext())
    assert WorldContext(geo_location=RefContext(name='x'))
    c1 = BioContext(cell_type=RefContext(name='x'))
    c2 = BioContext(cell_type=RefContext(name='y'))
    assert c1 != c2
    assert not c1 == c2
    assert c1 == c1
    assert not c1 != c1


def test_deprecated_cellular_location():
    stmt = Translocation(Agent('x'), 'HCN4 channel complex',
                         'pre-autophagosomal structure')
    assert stmt.from_location == 'HCN4 channel complex'
    assert stmt.to_location == 'pre-autophagosomal structure'
    stmt = Statement._from_json(stmt.to_json())
    assert stmt.from_location == 'HCN4 channel complex'
    assert stmt.to_location == 'pre-autophagosomal structure'


def test_stmt_type():
    stmt = Phosphorylation(None, Agent('x'))
    t = stmt_type(stmt)
    assert t == type(stmt), t
    mc = ModCondition('phosphorylation', 'S', '123')
    t = stmt_type(mc)
    assert t == 'ModCondition', t


def test_mk_str():
    stmt = Phosphorylation(None, Agent('x'))
    mk = stmt.matches_key()
    assert mk.startswith("(<class \'indra.statements.Phosphorylation\'>")


def test_ev_str():
    ev = Evidence(annotations={'a': 'b'}, text='test', source_api='test')
    str(ev)
    ev.__repr__()


@unittest.skip('Travis cannot draw the graph with pygraphviz.')
def test_draw_statements():
    stmt = Phosphorylation(None, Agent('x'))
    draw_stmt_graph([stmt])


def test_migration():
    m = Migration(Concept('migration'), 
                  QuantitativeState('person', 500, 'absolute'),
                  MovementContext(locations=[
                      {'location': RefContext('South Sudan'), 'role': 'origin'},
                      {'location': RefContext('Ethiopia'), 'role': 'destination'}],
                      time=TimeContext(text='beginning of 2019')))
    assert m
    assert isinstance(m, Migration)
    assert isinstance(m.delta, QuantitativeState)
    assert isinstance(m.context, MovementContext)
    assert m.context.locations[0]['location'].name == 'South Sudan'
    assert m.context.locations[0]['role'] == 'origin'
