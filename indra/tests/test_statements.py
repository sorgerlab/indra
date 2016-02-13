from indra.statements import *

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
    st1 = Phosphorylation(raf, mek, 'Phosphorylation', None)
    st2 = Phosphorylation(raf, mek, 'Phosphorylation', None)
    assert(st1.matches(st2))

def test_matches_key2():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek, 'Phosphorylation', None)
    st2 = Phosphorylation(raf, mek, 'Phosphorylation', None)
    assert(st1.matches_key() == st2.matches_key())

def test_not_matches():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek, 'Phosphorylation', None)
    st2 = Phosphorylation(raf, mek, 'PhosphorylationTyrosine', None)
    assert(not st1.matches(st2))

def test_not_matches_key():
    raf = Agent('Raf')
    mek = Agent('Mek')
    st1 = Phosphorylation(raf, mek, 'Phosphorylation', None)
    st2 = Phosphorylation(raf, mek, 'PhosphorylationTyrosine', None)
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
    st1 = Phosphorylation(src, nras1, 'PhosphorylationTyrosine', '32',
                          evidence=Evidence(text='foo'))
    st2 = Phosphorylation(src, nras2, 'Phosphorylation', None,
                          evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_selfmod():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = Autophosphorylation(nras1, 'PhosphorylationTyrosine', '32',
                          evidence=Evidence(text='foo'))
    st2 = Autophosphorylation(nras2, 'Phosphorylation', None,
                          evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_activityactivity():
    """Test matching of entities only, entities match only on name."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivityActivity(src, 'Kinase1', nras1, 'GtpBoundActivity1',
                           'increases1', evidence=Evidence(text='foo'))
    st2 = ActivityActivity(src, 'Kinase2', nras2, 'GtpBoundActivity2',
                           'increases2', evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_activitymod():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivityModification(nras1, 'PhosphorylationTyrosine', '32',
                               'increases1', 'GtpBoundActivity1',
                               evidence=Evidence(text='foo'))
    st1 = ActivityModification(nras2, 'Phosphorylation', None,
                               'increases2', 'GtpBoundActivity2',
                               evidence=Evidence(text='bar'))
    assert(st1.entities_match(st2))

def test_entities_match_activatingsub():
    """Test matching of entities only, entities match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = ActivatingSubstitution(nras1, 'G', '12', 'D', 'GtpBoundActivity1',
                                 'increases1', evidence=Evidence(text='foo'))
    st1 = ActivatingSubstitution(nras2, 'Q', '61', 'L', 'GtpBoundActivity2',
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

