from indra.preassembler import check_statements, check_sequence
from indra.preassembler import check_agent_mod
from indra.statements import Phosphorylation, Agent

def test_check_agent_mod_args():
    a = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a, ['PhosphorylationThreonine'], [185])
    assert(not failures)

def test_check_agent_mod_multiple():
    a = Agent('MAPK1', mods=['PhosphorylationThreonine', 'PhosphorylationTyrosine'],
              mod_sites=[185, 187], db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a)
    assert(not failures)

def test_check_agent_mod_wrong():
    a = Agent('MAPK1', mods=['PhosphorylationThreonine'],
              mod_sites=[186], db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a)
    assert(failures)

def test_check_agent_mod_wrong_multiple():
    a = Agent('MAPK1', mods=['PhosphorylationThreonine', 'PhosphorylationTyrosine'],
              mod_sites=[184, 186], db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a)
    assert(failures)

def test_check_agent_mod_wrong_multiple2():
    a = Agent('MAPK1', mods=['PhosphorylationThreonine', 'PhosphorylationTyrosine'],
              mod_sites=[185, 186], db_refs = {'UP': 'P28482'})
    failures = check_agent_mod(a)
    assert(failures)

def test_check_sequence_phos():
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1')
    stmt = Phosphorylation(e, s, 'PhosphorylationThreonine', 185)
    failures = check_sequence(stmt)
    assert(not failures)

def test_check_sequence_phos_wrong():
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1')
    stmt = Phosphorylation(e, s, 'PhosphorylationThreonine', 186)
    failures = check_sequence(stmt)
    assert(failures)

def test_check_sequence_phos_enz():
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1', mods=['PhosphorylationSerine'], mod_sites=[222], 
              db_refs = {'UP': 'Q02750'})
    stmt = Phosphorylation(e, s, 'PhosphorylationThreonine', 185)
    failures = check_sequence(stmt)
    assert(not failures)

def test_check_sequence_phos_enz_wrong():
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1', mods=['PhosphorylationSerine'], mod_sites=[221],
              db_refs = {'UP': 'Q02750'})
    stmt = Phosphorylation(e, s, 'PhosphorylationThreonine', 185)
    failures = check_sequence(stmt)
    assert(failures)

def test_check_statements():
    stmts = []
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1', mods=['PhosphorylationSerine'], mod_sites=[222], 
              db_refs = {'UP': 'Q02750'})
    stmt = Phosphorylation(e, s, 'PhosphorylationThreonine', 185)
    stmts.append(stmt)
    
    s = Agent('MAPK1', db_refs = {'UP': 'P28482'})
    e = Agent('MAP2K1', mods=['PhosphorylationSerine'], mod_sites=[221],
              db_refs = {'UP': 'Q02750'})
    stmt = Phosphorylation(e, s, 'PhosphorylationThreonine', 185)
    stmts.append(stmt)

    p, f = check_statements(stmts)
    assert(len(p) == 1)
    assert(len(f) == 1)
