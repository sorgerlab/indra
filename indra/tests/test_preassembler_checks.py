from indra.preassembler import check_statements, check_sequence
from indra.preassembler import check_agent_mod
from indra.statements import Phosphorylation, Agent, ModCondition

def test_check_agent_mod_args():
    a = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    mc = ModCondition('phosphorylation', 'threonine', '185')
    failures = check_agent_mod(a, [mc])
    assert(not failures)

def test_check_agent_mod_multiple():
    mc1 = ModCondition('phosphorylation', 'threonine', '185')
    mc2 = ModCondition('phosphorylation', 'tyrosine', '187')
    a = Agent('MAPK1', mods=[mc1, mc2], db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a)
    assert(not failures)

def test_check_agent_mod_wrong():
    mc = ModCondition('phosphorylation', 'threonine', '186')
    a = Agent('MAPK1', mods=[mc], db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a)
    assert(failures)

def test_check_agent_mod_wrong_multiple():
    mc1 = ModCondition('phosphorylation', 'threonine', '184')
    mc2 = ModCondition('phosphorylation', 'tyrosine', '186')
    a = Agent('MAPK1', mods=[mc1, mc2], db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a)
    assert(failures)

def test_check_agent_mod_wrong_multiple2():
    mc1 = ModCondition('phosphorylation', 'threonine', '185')
    mc2 = ModCondition('phosphorylation', 'tyrosine', '186')
    a = Agent('MAPK1', mods=[mc1, mc2], db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a)
    assert(failures)

def test_check_sequence_phos():
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1')
    mc1 = ModCondition('phosphorylation', 'threonine', '185')
    stmt = Phosphorylation(e, s, mc1)
    failures = check_sequence(stmt)
    assert(not failures)

def test_check_sequence_phos_wrong():
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1')
    mc1 = ModCondition('phosphorylation', 'threonine', '186')
    stmt = Phosphorylation(e, s, mc1)
    failures = check_sequence(stmt)
    assert(failures)

def test_check_sequence_phos_enz():
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    mc1 = ModCondition('phosphorylation', 'serine', '222')
    mc2 = ModCondition('phosphorylation', 'threoninine', '185')
    e = Agent('MAP2K1', mc1, db_refs = {'UP': 'Q02750'})
    stmt = Phosphorylation(e, s, mc2)
    failures = check_sequence(stmt)
    assert(not failures)

def test_check_sequence_phos_enz_wrong():
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    mc1 = ModCondition('phosphorylation', 'serine', '221')
    mc2 = ModCondition('phosphorylation', 'threoninine', '185')
    e = Agent('MAP2K1', mods=[mc1],
              db_refs = {'UP': 'Q02750'})
    stmt = Phosphorylation(e, s, mc2)
    failures = check_sequence(stmt)
    assert(failures)

def test_check_statements():
    stmts = []
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    mc1 = ModCondition('phosphorylation', 'serine', '222')
    mc2 = ModCondition('phosphorylation', 'threoninine', '185')
    e = Agent('MAP2K1', mods=[mc1],
              db_refs = {'UP': 'Q02750'})
    stmt = Phosphorylation(e, s, mc2)
    stmts.append(stmt)

    mc1 = ModCondition('phosphorylation', 'serine', '221')
    mc2 = ModCondition('phosphorylation', 'threoninine', '185')
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1', mods=[mc1], db_refs = {'UP': 'Q02750'})
    stmt = Phosphorylation(e, s, mc2)
    stmts.append(stmt)

    p, f = check_statements(stmts)
    assert(len(p) == 1)
    assert(len(f) == 1)
