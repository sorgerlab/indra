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

def test_agent_name_match():
    """Agents match only on name."""
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    assert(nras1.entity_matches(nras2))

def test_agents_match_phos():
    """Test matching of agents only, agents match only on name."""
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    nras1 = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nras2 = Agent('NRAS', db_refs = {'HGNC': 'dummy'})
    st1 = Phosphorylation(src, nras1, 'PhosphorylationTyrosine', '32',
                          evidence=Evidence(text='foo'))
    st2 = Phosphorylation(src, nras2, 'Phosphorylation', None,
                          evidence=Evidence(text='bar'))
    assert(st1.agents_match(st2))

