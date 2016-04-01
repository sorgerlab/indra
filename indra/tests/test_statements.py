from indra.statements import *
import pkg_resources
import os
from indra.preassembler.hierarchy_manager import HierarchyManager
from nose.tools import raises

# Load the hierarchy manager data
ent_path = os.path.join('preassembler', 'entity_hierarchy.rdf')
ent_file = pkg_resources.resource_filename('indra', ent_path)
eh = HierarchyManager(ent_file)

mod_path = os.path.join('preassembler', 'modification_hierarchy.rdf')
mod_file = pkg_resources.resource_filename('indra', mod_path)
mh = HierarchyManager(mod_file)

# Argument checking for ActivityModifications ----------------------------

def test_activitymod_sitelist_of_ints():
    """Check that mod positions specified as ints get promoted to strings."""
    st = ActivityModification(Agent('MAP2K1'),
                              [ModCondition('phosphorylation', 'serine', 218),
                               ModCondition('phosphorylation', 'serine', 222)],
                              'increases', 'KinaseActivity')
    assert isinstance(st.mod[0].position, basestring)
    assert isinstance(st.mod[1].position, basestring)

def testactivitymod_string_string():
    """Check that string mod position is preserved"""
    st = ActivityModification(Agent('MAP2K1'),
                              [ModCondition('phosphorylation', 'serine', '218')],
                              'increases', 'KinaseActivity')
    assert isinstance(st.mod[0].position, basestring)

def testactivitymod_string_none():
    """Check that None mod position is preserved"""
    st = ActivityModification(Agent('MAP2K1'),
                              [ModCondition('phosphorylation', 'serine', None)],
                              'increases', 'KinaseActivity')
    assert (st.mod[0].position is None)

def testactivitymod_nolist():
    """Make sure mod is correctly turned into a list if it's 
    a single ModCondition"""
    mc = ModCondition('phosphorylation')
    st = ActivityModification(Agent('MAP2K1'), mc,
                              'increases', 'KinaseActivity')
    assert isinstance(st.mod, list)

# Checking for exact matching (except Evidence) between Agents/stmts ---------

def test_matches():
    ras = Agent('Ras')
    raf = Agent('Raf')
    st1 = Complex([ras, raf])
    st2 = Complex([ras, raf])
    assert(st1.matches(st2))

def test_matches_key():
    ras = Agent('Ras')
    raf = Agent('Raf')
    st1 = Complex([ras, raf])
    st2 = Complex([ras, raf])
    assert(st1.matches_key() == st2.matches_key())

def test_matches2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek)
    assert(st1.matches(st2))

def test_matches_key2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek)
    assert(st1.matches_key() == st2.matches_key())

def test_not_matches():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek, 'tyrosine')
    assert(not st1.matches(st2))

def test_not_matches_key():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek)
    st2 = Phosphorylation(raf, mek, 'tyrosine')
    assert(st1.matches_key() != st2.matches_key())

def test_matches_dbrefs():
    hras1 = Agent('HRAS', db_refs={'hgnc': 1111})
    hras2 = Agent('HRAS', db_refs={'hgnc': 9999})
    assert(hras1.matches(hras2))

def test_matches_key_dbrefs():
    hras1 = Agent('HRAS', db_refs={'hgnc': 1111})
    hras2 = Agent('HRAS', db_refs={'hgnc': 9999})
    assert(hras1.matches_key() == hras2.matches_key())

def test_matches_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    assert(hras1.matches(hras2))

def test_matches_key_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    assert(hras1.matches_key() == hras2.matches_key())

def test_not_matches_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('RAF1'), True)])
    assert(not hras1.matches(hras2))

def test_not_matches_key_bound():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches_key() != hras2.matches_key())

def test_not_matches_bound2():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), False)])
    assert(not hras1.matches(hras2))

def test_not_matches_key_bound2():
    hras1 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS',
        bound_conditions=[BoundCondition(Agent('BRAF'), False)])
    assert(hras1.matches_key() != hras2.matches_key())

def test_matches_bound_multiple():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches(hras2))

def test_matches_key_bound_multiple():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches_key() == hras2.matches_key())

def test_matches_bound_multiple_order():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('RAF1'), True),
                                            BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches(hras2))

def test_matches_key_bound_multiple_order():
    hras1 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('RAF1'), True),
                                            BoundCondition(Agent('BRAF'), True)])
    hras2 = Agent('HRAS', bound_conditions=[BoundCondition(Agent('BRAF'), True),
                                            BoundCondition(Agent('RAF1'), True)])
    assert(hras1.matches_key() == hras2.matches_key())

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

def test_matches_activityactivity():
    """Test matching of entities only, entities match only on name."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivityActivity(src, 'Kinase1', 'increases1',
                           nras1, 'GtpBoundActivity1',
                           evidence=Evidence(text='foo'))
    st2 = ActivityActivity(src, 'Kinase1', 'increases1',
                           nras1, 'GtpBoundActivity1',
                           evidence=Evidence(text='bar'))
    st3 = ActivityActivity(src, 'Kinase2', 'increases2',
                           nras2, 'GtpBoundActivity2', 
                           evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))

def test_matches_activitymod():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivityModification(nras1, ModCondition('phosphorylation',
                                                   'tyrosine', '32'),
                               'increases1', 'GtpBoundActivity1',
                               evidence=Evidence(text='foo'))
    st2 = ActivityModification(nras1, ModCondition('phosphorylation',
                                                   'tyrosine', '32'),
                               'increases1', 'GtpBoundActivity1',
                               evidence=Evidence(text='bar'))
    st3 = ActivityModification(nras2, ModCondition('phosphorylation'),
                               'increases2', 'GtpBoundActivity2',
                               evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))

def test_matches_activatingsub():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivatingSubstitution(nras1, 'G', '12', 'D', 'GtpBoundActivity1',
                                 'increases1', evidence=Evidence(text='foo'))
    st2 = ActivatingSubstitution(nras1, 'G', '12', 'D', 'GtpBoundActivity1',
                                 'increases1', evidence=Evidence(text='bar'))
    st3 = ActivatingSubstitution(nras2, 'Q', '61', 'L', 'GtpBoundActivity2',
                                 'increases2', evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))

def test_matches_rasgef():
    """Test matching of entities only, entities match only on name."""
    sos1 = Agent('SOS1', db_refs = {'HGNC': 'sos1'})
    sos2 = Agent('SOS1', db_refs = {'HGNC': 'sos2'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = RasGef(sos1, 'GtpBoundActivity1', nras1,
                 evidence=Evidence(text='foo'))
    st2 = RasGef(sos1, 'GtpBoundActivity1', nras1,
                 evidence=Evidence(text='bar'))
    st3 = RasGef(sos2, 'GtpBoundActivity2', nras2,
                 evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))

def test_matches_rasgap():
    rasa1 = Agent('RASA1', db_refs = {'HGNC': 'rasa1'})
    rasa2 = Agent('RASA1', db_refs = {'HGNC': 'rasa2'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = RasGap(rasa1, 'GtpBoundActivity1', nras1,
                 evidence=Evidence(text='foo'))
    st2 = RasGap(rasa1, 'GtpBoundActivity1', nras1,
                 evidence=Evidence(text='bar'))
    st3 = RasGap(rasa2, 'GtpBoundActivity2', nras2,
                 evidence=Evidence(text='bar'))
    assert(st1.matches(st2))
    assert(not st1.matches(st3))

def test_matches_complex():
    ksr1 = Agent('KSR1', db_refs = {'HGNC': 'ksr1'})
    ksr2 = Agent('KSR1', db_refs = {'HGNC': 'ksr2'})
    braf1 = Agent('BRAF', db_refs = {'HGNC': 'braf1'})
    braf2 = Agent('BRAF', db_refs = {'HGNC': 'braf2'})
    map2k1 = Agent('MAP2K1', db_refs = {'HGNC': 'map2k1'})
    map2k2 = Agent('MAP2K1', db_refs = {'HGNC': 'map2k2'})
    st1 = Complex([ksr1, braf1, map2k1], evidence=Evidence(text='foo'))
    st2 = Complex([ksr1, braf1, map2k1], evidence=Evidence(text='bar'))
    assert(st1.matches(st2))


# Entity matching between statements ----------------------------------------
def test_agent_entity_match():
    """Agents match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    assert(nras1.entity_matches(nras2))

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

def test_entities_match_selfmod():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = Autophosphorylation(nras1, 'tyrosine', '32',
                          evidence=Evidence(text='foo'))
    st2 = Autophosphorylation(nras2,
                          evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_activityactivity():
    """Test matching of entities only, entities match only on name."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivityActivity(src, 'Kinase1', 'increases1',
                           nras1, 'GtpBoundActivity1',
                           evidence=Evidence(text='foo'))
    st2 = ActivityActivity(src, 'Kinase2', 'increases2',
                           nras2, 'GtpBoundActivity2', 
                           evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_activitymod():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivityModification(nras1, ModCondition('phosphorylation',
                                                   'tyrosine', '32'),
                               'increases1', 'GtpBoundActivity1',
                               evidence=Evidence(text='foo'))
    st2 = ActivityModification(nras2, ModCondition('phosphorylation'),
                               'increases2', 'GtpBoundActivity2',
                               evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_activatingsub():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivatingSubstitution(nras1, 'G', '12', 'D', 'GtpBoundActivity1',
                                 'increases1', evidence=Evidence(text='foo'))
    st2 = ActivatingSubstitution(nras2, 'Q', '61', 'L', 'GtpBoundActivity2',
                                 'increases2', evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_rasgef():
    """Test matching of entities only, entities match only on name."""
    sos1 = Agent('SOS1', db_refs = {'HGNC': 'sos1'})
    sos2 = Agent('SOS1', db_refs = {'HGNC': 'sos2'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = RasGef(sos1, 'GtpBoundActivity1', nras1,
                 evidence=Evidence(text='foo'))
    st2 = RasGef(sos2, 'GtpBoundActivity2', nras2,
                 evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_rasgap():
    """Test matching of entities only, entities match only on name."""
    rasa1 = Agent('RASA1', db_refs = {'HGNC': 'rasa1'})
    rasa2 = Agent('RASA1', db_refs = {'HGNC': 'rasa2'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = RasGap(rasa1, 'GtpBoundActivity1', nras1,
                 evidence=Evidence(text='foo'))
    st2 = RasGap(rasa2, 'GtpBoundActivity2', nras2,
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
    assert(st1.entities_match(st2))

def test_agent_superfamily_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    ras = Agent('RAS', db_refs = {'FA': '03663'})
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    assert nras.refinement_of(ras, eh, mh)
    assert not ras.refinement_of(nras, eh, mh)
    # The top-level list should contain only one statement, the gene-level
    # one, supported by the family one.

def test_agent_boundcondition_refinement():
    """A gene-level statement should be supported by a family-level
    statement."""
    bc1 = BoundCondition(Agent('BRAF', db_refs = {'HGNC': 'braf1'}), True)
    bc2 = BoundCondition(Agent('RAF1', db_refs = {'HGNC': 'braf2'}), True)
    bc3 = BoundCondition(Agent('RAF1', db_refs = {'HGNC': 'braf2'}), False)
    bc4 = BoundCondition(Agent('RAF', db_refs = {'HGNC': 'raf'}), True)

    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'}, bound_conditions=[bc1])
    nras2 = Agent('NRAS', db_refs = {'HGNC': '7989'}, bound_conditions=[bc2])
    nras3 = Agent('NRAS', db_refs = {'HGNC': '7989'}, bound_conditions=[bc3])
    nras4 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras5 = Agent('NRAS', db_refs = {'HGNC': '7989'},
                  bound_conditions=[bc4])

    # nras1 (bound to BRAF)
    assert not nras2.refinement_of(nras1, eh, mh)
    assert not nras3.refinement_of(nras1, eh, mh)
    assert not nras4.refinement_of(nras1, eh, mh)
    assert not nras5.refinement_of(nras1, eh, mh)
    # nras2 (bound to CRAF)
    assert not nras1.refinement_of(nras2, eh, mh)
    assert not nras3.refinement_of(nras2, eh, mh) # Not bound condition
    assert not nras4.refinement_of(nras2, eh, mh)
    assert not nras5.refinement_of(nras2, eh, mh)
    # nras3 (not bound to CRAF)
    assert not nras1.refinement_of(nras3, eh, mh)
    assert not nras2.refinement_of(nras3, eh, mh) # Not bound condition
    assert not nras4.refinement_of(nras3, eh, mh)
    assert not nras5.refinement_of(nras3, eh, mh)
    # nras4 (no bound condition)
    assert nras1.refinement_of(nras4, eh, mh)
    assert nras2.refinement_of(nras4, eh, mh)
    assert nras3.refinement_of(nras4, eh, mh)
    assert nras5.refinement_of(nras4, eh, mh)
    # nras5 (RAF family bound condition)
    assert nras1.refinement_of(nras5, eh, mh)
    assert nras2.refinement_of(nras5, eh, mh)
    assert not nras3.refinement_of(nras5, eh, mh)
    assert not nras4.refinement_of(nras5, eh, mh)

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
    assert mek2.refinement_of(mek1, eh, mh)
    assert mek3.refinement_of(mek1, eh, mh)
    assert mek4.refinement_of(mek1, eh, mh)
    assert mek5.refinement_of(mek1, eh, mh)
    assert mek6.refinement_of(mek1, eh, mh)
    assert mek7.refinement_of(mek1, eh, mh)
    assert mek8.refinement_of(mek1, eh, mh)
    # mek2
    assert not mek1.refinement_of(mek2, eh, mh)
    assert not mek3.refinement_of(mek2, eh, mh) # Different site
    assert mek4.refinement_of(mek2, eh, mh)
    assert not mek5.refinement_of(mek2, eh, mh) # Cross-relationship
    assert mek6.refinement_of(mek2, eh, mh)
    assert not mek7.refinement_of(mek2, eh, mh) # Different site
    assert mek8.refinement_of(mek2, eh, mh)
    # mek3
    assert not mek1.refinement_of(mek3, eh, mh)
    assert not mek2.refinement_of(mek3, eh, mh)
    assert mek4.refinement_of(mek3, eh, mh)
    assert not mek5.refinement_of(mek3, eh, mh)
    assert not mek6.refinement_of(mek3, eh, mh)
    assert mek7.refinement_of(mek3, eh, mh)
    assert mek8.refinement_of(mek3, eh, mh)
    # mek4
    assert not mek1.refinement_of(mek4, eh, mh)
    assert not mek2.refinement_of(mek4, eh, mh)
    assert not mek3.refinement_of(mek4, eh, mh)
    assert not mek5.refinement_of(mek4, eh, mh)
    assert not mek6.refinement_of(mek4, eh, mh)
    assert not mek7.refinement_of(mek4, eh, mh)
    assert mek8.refinement_of(mek4, eh, mh)
    # mek5
    assert not mek1.refinement_of(mek5, eh, mh)
    assert not mek2.refinement_of(mek5, eh, mh)
    assert not mek3.refinement_of(mek5, eh, mh)
    assert not mek4.refinement_of(mek5, eh, mh)
    assert mek6.refinement_of(mek5, eh, mh)
    assert mek7.refinement_of(mek5, eh, mh)
    assert mek8.refinement_of(mek5, eh, mh)
    # mek6
    assert not mek1.refinement_of(mek6, eh, mh)
    assert not mek2.refinement_of(mek6, eh, mh)
    assert not mek3.refinement_of(mek6, eh, mh)
    assert not mek4.refinement_of(mek6, eh, mh)
    assert not mek5.refinement_of(mek6, eh, mh)
    assert not mek7.refinement_of(mek6, eh, mh)
    assert mek8.refinement_of(mek6, eh, mh)
    # mek7
    assert not mek1.refinement_of(mek7, eh, mh)
    assert not mek2.refinement_of(mek7, eh, mh)
    assert not mek3.refinement_of(mek7, eh, mh)
    assert not mek4.refinement_of(mek7, eh, mh)
    assert not mek5.refinement_of(mek7, eh, mh)
    assert not mek6.refinement_of(mek7, eh, mh)
    assert mek8.refinement_of(mek7, eh, mh)
    # mek8
    assert not mek1.refinement_of(mek8, eh, mh)
    assert not mek2.refinement_of(mek8, eh, mh)
    assert not mek3.refinement_of(mek8, eh, mh)
    assert not mek4.refinement_of(mek8, eh, mh)
    assert not mek5.refinement_of(mek8, eh, mh)
    assert not mek6.refinement_of(mek8, eh, mh)
    assert not mek7.refinement_of(mek8, eh, mh)

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
    assert p2.refinement_of(p1, eh, mh)
    assert p3.refinement_of(p1, eh, mh)
    assert p4.refinement_of(p1, eh, mh)
    assert p5.refinement_of(p1, eh, mh)
    assert p6.refinement_of(p1, eh, mh)
    # p2
    assert not p1.refinement_of(p2, eh, mh)
    assert not p3.refinement_of(p2, eh, mh)
    assert not p4.refinement_of(p2, eh, mh)
    assert p5.refinement_of(p2, eh, mh)
    assert not p6.refinement_of(p2, eh, mh)
    # p3
    assert not p1.refinement_of(p3, eh, mh)
    assert not p2.refinement_of(p3, eh, mh)
    assert not p4.refinement_of(p3, eh, mh)
    assert not p5.refinement_of(p3, eh, mh)
    assert p6.refinement_of(p3, eh, mh)
    # p4
    assert not p1.refinement_of(p4, eh, mh)
    assert not p2.refinement_of(p4, eh, mh)
    assert not p3.refinement_of(p4, eh, mh)
    assert p5.refinement_of(p4, eh, mh)
    assert p6.refinement_of(p4, eh, mh)
    # p5
    assert not p1.refinement_of(p5, eh, mh)
    assert not p2.refinement_of(p5, eh, mh)
    assert not p3.refinement_of(p5, eh, mh)
    assert not p4.refinement_of(p5, eh, mh)
    assert not p6.refinement_of(p5, eh, mh)
    # p6
    assert not p1.refinement_of(p6, eh, mh)
    assert not p2.refinement_of(p6, eh, mh)
    assert not p3.refinement_of(p6, eh, mh)
    assert not p4.refinement_of(p6, eh, mh)
    assert not p5.refinement_of(p6, eh, mh)

def test_autophosphorylation_modification_refinement():
    braf = Agent('BRAF', db_refs = {'HGNC': 'braf'})
    p1 = Autophosphorylation(braf,)
    p2 = Autophosphorylation(braf, position='218')
    p3 = Autophosphorylation(braf, position='222')
    p4 = Autophosphorylation(braf, 'serine')
    p5 = Autophosphorylation(braf, 'serine', '218')
    p6 = Autophosphorylation(braf, 'serine', '222')

    # p1
    assert p2.refinement_of(p1, eh, mh)
    assert p3.refinement_of(p1, eh, mh)
    assert p4.refinement_of(p1, eh, mh)
    assert p5.refinement_of(p1, eh, mh)
    assert p6.refinement_of(p1, eh, mh)
    # p2
    assert not p1.refinement_of(p2, eh, mh)
    assert not p3.refinement_of(p2, eh, mh)
    assert not p4.refinement_of(p2, eh, mh)
    assert p5.refinement_of(p2, eh, mh)
    assert not p6.refinement_of(p2, eh, mh)
    # p3
    assert not p1.refinement_of(p3, eh, mh)
    assert not p2.refinement_of(p3, eh, mh)
    assert not p4.refinement_of(p3, eh, mh)
    assert not p5.refinement_of(p3, eh, mh)
    assert p6.refinement_of(p3, eh, mh)
    # p4
    assert not p1.refinement_of(p4, eh, mh)
    assert not p2.refinement_of(p4, eh, mh)
    assert not p3.refinement_of(p4, eh, mh)
    assert p5.refinement_of(p4, eh, mh)
    assert p6.refinement_of(p4, eh, mh)
    # p5
    assert not p1.refinement_of(p5, eh, mh)
    assert not p2.refinement_of(p5, eh, mh)
    assert not p3.refinement_of(p5, eh, mh)
    assert not p4.refinement_of(p5, eh, mh)
    assert not p6.refinement_of(p5, eh, mh)
    # p6
    assert not p1.refinement_of(p6, eh, mh)
    assert not p2.refinement_of(p6, eh, mh)
    assert not p3.refinement_of(p6, eh, mh)
    assert not p4.refinement_of(p6, eh, mh)
    assert not p5.refinement_of(p6, eh, mh)

def test_activityactivity_modification_refinement():
    raf = Agent('RAF')
    braf = Agent('BRAF')
    mek = Agent('MEK')
    mek1 = Agent('MAP2K1')

    st1 = ActivityActivity(raf, 'KinaseActivity', 'increases',
                           mek, 'KinaseActivity')
    st2 = ActivityActivity(braf, 'KinaseActivity', 'increases',
                           mek, 'KinaseActivity')
    st3 = ActivityActivity(raf, 'KinaseActivity', 'increases',
                           mek1, 'KinaseActivity')
    st4 = ActivityActivity(braf, 'KinaseActivity', 'increases',
                           mek1, 'KinaseActivity')
    st5 = ActivityActivity(braf, 'KinaseActivity', 'increasesX',
                           mek1, 'KinaseActivity')
    # st1
    assert st2.refinement_of(st1, eh, mh)
    assert st3.refinement_of(st1, eh, mh)
    assert st4.refinement_of(st1, eh, mh)
    assert not st5.refinement_of(st1, eh, mh)
    # st2
    assert not st1.refinement_of(st2, eh, mh)
    assert not st3.refinement_of(st2, eh, mh)
    assert st4.refinement_of(st2, eh, mh)
    assert not st5.refinement_of(st2, eh, mh)
    # st3
    assert not st1.refinement_of(st3, eh, mh)
    assert not st2.refinement_of(st3, eh, mh)
    assert st4.refinement_of(st3, eh, mh)
    assert not st5.refinement_of(st3, eh, mh)
    # st4
    assert not st1.refinement_of(st4, eh, mh)
    assert not st2.refinement_of(st4, eh, mh)
    assert not st3.refinement_of(st4, eh, mh)
    assert not st5.refinement_of(st4, eh, mh)
    # st5
    assert not st1.refinement_of(st5, eh, mh)
    assert not st2.refinement_of(st5, eh, mh)
    assert not st3.refinement_of(st5, eh, mh)
    assert not st4.refinement_of(st5, eh, mh)

def test_activitymod_refinement():
    mek_fam = Agent('MEK')
    mek1 = Agent('MAP2K1')
    p1 = ActivityModification(mek_fam, ModCondition('phosphorylation'),
                              'increases', 'KinaseActivity')
    p2 = ActivityModification(mek_fam, ModCondition('phosphorylation',
                                                    'serine', '218'),
                              'increases', 'KinaseActivity')
    p3 = ActivityModification(mek1, ModCondition('phosphorylation'),
                              'increases', 'KinaseActivity')
    p4 = ActivityModification(mek1, ModCondition('phosphorylation', 'serine'),
                              'increases', 'KinaseActivity')
    p5 = ActivityModification(mek1, ModCondition('phosphorylation',
                                                 'serine', '218'),
                              'increases', 'KinaseActivity')
    p6 = ActivityModification(mek1, ModCondition('phosphorylation',
                                                 'serine', '222'),
                              'increases', 'KinaseActivity')
    p7 = ActivityModification(mek1,
                            [ModCondition('phosphorylation', 'serine', '218'),
                             ModCondition('phosphorylation', 'serine', '222')],
                            'increases', 'KinaseActivity')
    # p1
    assert p2.refinement_of(p1, eh, mh)
    assert p3.refinement_of(p1, eh, mh)
    assert p4.refinement_of(p1, eh, mh)
    assert p5.refinement_of(p1, eh, mh)
    assert p6.refinement_of(p1, eh, mh)
    assert p7.refinement_of(p1, eh, mh)
    # p2
    assert not p1.refinement_of(p2, eh, mh)
    assert not p3.refinement_of(p2, eh, mh)
    assert not p4.refinement_of(p2, eh, mh)
    assert p5.refinement_of(p2, eh, mh)
    assert not p6.refinement_of(p2, eh, mh)
    assert p7.refinement_of(p2, eh, mh)
    # p3
    assert not p1.refinement_of(p3, eh, mh)
    assert not p2.refinement_of(p3, eh, mh)
    assert p4.refinement_of(p3, eh, mh)
    assert p5.refinement_of(p3, eh, mh)
    assert p6.refinement_of(p3, eh, mh)
    assert p7.refinement_of(p3, eh, mh)
    # p4
    assert not p1.refinement_of(p4, eh, mh)
    assert not p2.refinement_of(p4, eh, mh)
    assert not p3.refinement_of(p4, eh, mh)
    assert p5.refinement_of(p4, eh, mh)
    assert p6.refinement_of(p4, eh, mh)
    assert p7.refinement_of(p4, eh, mh)
    # p5
    assert not p1.refinement_of(p5, eh, mh)
    assert not p2.refinement_of(p5, eh, mh)
    assert not p3.refinement_of(p5, eh, mh)
    assert not p4.refinement_of(p5, eh, mh)
    assert not p6.refinement_of(p5, eh, mh)
    assert p7.refinement_of(p5, eh, mh)
    # p6
    assert not p1.refinement_of(p6, eh, mh)
    assert not p2.refinement_of(p6, eh, mh)
    assert not p3.refinement_of(p6, eh, mh)
    assert not p4.refinement_of(p6, eh, mh)
    assert not p5.refinement_of(p6, eh, mh)
    assert p7.refinement_of(p6, eh, mh)
    # p7
    assert not p1.refinement_of(p7, eh, mh)
    assert not p2.refinement_of(p7, eh, mh)
    assert not p3.refinement_of(p7, eh, mh)
    assert not p4.refinement_of(p7, eh, mh)
    assert not p5.refinement_of(p7, eh, mh)
    assert not p6.refinement_of(p7, eh, mh)

def test_activatingsub_family_refinement():
    ras = Agent('RAS')
    kras = Agent('KRAS')
    nras = Agent('NRAS')
    st1 = ActivatingSubstitution(ras, 'G', '12', 'D', 'activity', 'increases')
    st2 = ActivatingSubstitution(kras, 'G', '12', 'D', 'activity', 'increases')
    st3 = ActivatingSubstitution(nras, 'G', '12', 'D', 'activity', 'increases')
    st4 = ActivatingSubstitution(kras, 'G', '12', 'D', 'activity', 'increasesX')
    # st1
    assert st2.refinement_of(st1, eh, mh)
    assert st3.refinement_of(st1, eh, mh)
    assert not st4.refinement_of(st1, eh, mh)
    # st2
    assert not st1.refinement_of(st2, eh, mh)
    assert not st3.refinement_of(st2, eh, mh)
    assert not st4.refinement_of(st2, eh, mh)
    # st3
    assert not st1.refinement_of(st3, eh, mh)
    assert not st2.refinement_of(st3, eh, mh)
    assert not st4.refinement_of(st3, eh, mh)
    # st4
    assert not st1.refinement_of(st4, eh, mh)
    assert not st2.refinement_of(st4, eh, mh)
    assert not st3.refinement_of(st4, eh, mh)

def test_rasgef_family_refinement():
    sos = Agent('SOS')
    sos1 = Agent('SOS1')
    ras = Agent('RAS')
    kras = Agent('KRAS')
    # Statements
    st1 = RasGef(sos, 'activity', ras)
    st2 = RasGef(sos1, 'activity', ras)
    st3 = RasGef(sos, 'activity', kras)
    st4 = RasGef(sos1, 'activity', kras)
    st5 = RasGef(sos1, 'different_activity', kras)
    # st1
    assert st2.refinement_of(st1, eh, mh)
    assert st3.refinement_of(st1, eh, mh)
    assert st4.refinement_of(st1, eh, mh)
    assert not st5.refinement_of(st1, eh, mh)
    # st2
    assert not st1.refinement_of(st2, eh, mh)
    assert not st3.refinement_of(st2, eh, mh)
    assert st4.refinement_of(st2, eh, mh)
    assert not st5.refinement_of(st2, eh, mh)
    # st3
    assert not st1.refinement_of(st3, eh, mh)
    assert not st2.refinement_of(st3, eh, mh)
    assert st4.refinement_of(st3, eh, mh)
    assert not st5.refinement_of(st3, eh, mh)
    # st4
    assert not st1.refinement_of(st4, eh, mh)
    assert not st2.refinement_of(st4, eh, mh)
    assert not st3.refinement_of(st4, eh, mh)
    assert not st5.refinement_of(st4, eh, mh)
    # st5
    assert not st1.refinement_of(st5, eh, mh)
    assert not st2.refinement_of(st5, eh, mh)
    assert not st3.refinement_of(st5, eh, mh)
    assert not st4.refinement_of(st5, eh, mh)

def test_rasgap_family_refinement():
    rasa = Agent('RASA')
    rasa1 = Agent('RASA1')
    ras = Agent('RAS')
    kras = Agent('KRAS')
    # Statements
    st1 = RasGap(rasa, 'activity', ras)
    st2 = RasGap(rasa1, 'activity', ras)
    st3 = RasGap(rasa, 'activity', kras)
    st4 = RasGap(rasa1, 'activity', kras)
    st5 = RasGap(rasa1, 'different_activity', kras)
    # st1
    assert st2.refinement_of(st1, eh, mh)
    assert st3.refinement_of(st1, eh, mh)
    assert st4.refinement_of(st1, eh, mh)
    assert not st5.refinement_of(st1, eh, mh)
    # st2
    assert not st1.refinement_of(st2, eh, mh)
    assert not st3.refinement_of(st2, eh, mh)
    assert st4.refinement_of(st2, eh, mh)
    assert not st5.refinement_of(st2, eh, mh)
    # st3
    assert not st1.refinement_of(st3, eh, mh)
    assert not st2.refinement_of(st3, eh, mh)
    assert st4.refinement_of(st3, eh, mh)
    assert not st5.refinement_of(st3, eh, mh)
    # st4
    assert not st1.refinement_of(st4, eh, mh)
    assert not st2.refinement_of(st4, eh, mh)
    assert not st3.refinement_of(st4, eh, mh)
    assert not st5.refinement_of(st4, eh, mh)
    # st5
    assert not st1.refinement_of(st5, eh, mh)
    assert not st2.refinement_of(st5, eh, mh)
    assert not st3.refinement_of(st5, eh, mh)
    assert not st4.refinement_of(st5, eh, mh)

def test_complex_family_refinement():
    raf = Agent('RAF')
    braf = Agent('BRAF')
    raf1 = Agent('RAF1')
    mek = Agent('MEK')
    mek1 = Agent('MAP2K1')

    st1 = Complex([raf, mek])
    st2 = Complex([braf, mek])
    st3 = Complex([raf, mek1])
    st4 = Complex([braf, mek1])
    st5 = Complex([braf, raf1])

    # st1
    assert st2.refinement_of(st1, eh, mh)
    assert st3.refinement_of(st1, eh, mh)
    assert st4.refinement_of(st1, eh, mh)
    assert not st5.refinement_of(st1, eh, mh)
    # st2
    assert not st1.refinement_of(st2, eh, mh)
    assert not st3.refinement_of(st2, eh, mh)
    assert st4.refinement_of(st2, eh, mh)
    assert not st5.refinement_of(st2, eh, mh)
    # st3
    assert not st1.refinement_of(st3, eh, mh)
    assert not st2.refinement_of(st3, eh, mh)
    assert st4.refinement_of(st3, eh, mh)
    assert not st5.refinement_of(st3, eh, mh)
    # st4
    assert not st1.refinement_of(st4, eh, mh)
    assert not st2.refinement_of(st4, eh, mh)
    assert not st3.refinement_of(st4, eh, mh)
    assert not st5.refinement_of(st4, eh, mh)
    # st5
    assert not st1.refinement_of(st5, eh, mh)
    assert not st2.refinement_of(st5, eh, mh)
    assert not st3.refinement_of(st5, eh, mh)
    assert not st4.refinement_of(st5, eh, mh)

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

def test_valid_residue():
    assert(get_valid_residue('serine') == 'S')
    assert(get_valid_residue('ser') == 'S')
    assert(get_valid_residue('Serine') == 'S')
    assert(get_valid_residue('SERINE') == 'S')

# TODO expand tests to also check for things that should NOT match (different
# agent names)

# TODO Expand match tests to 

