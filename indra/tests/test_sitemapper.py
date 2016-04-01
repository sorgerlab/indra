from indra.preassembler.sitemapper import default_mapper as sm, MappedStatement
from indra.statements import *

def test_check_agent_mod():
    mapk1_valid = Agent('MAPK1',
                        mods=[ModCondition('phosphorylation', 'T', '185'),
                              ModCondition('phosphorylation', 'Y', '187')],
                        db_refs={'UP': 'P28482'})
    res_valid = sm.map_agent_sites(mapk1_valid)
    assert len(res_valid) == 2
    assert res_valid[0] == []
    assert isinstance(res_valid[1], Agent)
    assert res_valid[1].matches(mapk1_valid)

    mapk1_invalid = Agent('MAPK1',
                          mods=[ModCondition('phosphorylation', 'T', '183'),
                                ModCondition('phosphorylation', 'Y', '185')],
                          db_refs={'UP': 'P28482'})
    res_invalid = sm.map_agent_sites(mapk1_invalid)
    assert len(res_invalid) == 2
    assert isinstance(res_invalid[0], list)
    assert isinstance(res_invalid[1], Agent)
    invalid_sites = res_invalid[0]
    assert len(invalid_sites) == 2
    map183 = invalid_sites[0]
    assert map183[0] == ('MAPK1', 'T', '183')
    assert len(map183[1]) == 3
    assert map183[1][0] == 'T'
    assert map183[1][1] == '185'
    map185 = invalid_sites[1]
    assert map185[0] == ('MAPK1', 'Y', '185')
    assert len(map185[1]) == 3
    assert map185[1][0] == 'Y'
    assert map185[1][1] == '187'
    new_agent = res_invalid[1]
    assert len(new_agent.mods) == 2
    assert new_agent.mods[0].matches(ModCondition('phosphorylation', 'T', '185'))
    assert new_agent.mods[1].matches(ModCondition('phosphorylation', 'Y', '187'))

    # TODO Test a site that is invalid but not found in the site map

def test_site_map_complex():
    mapk1_invalid = Agent('MAPK1',
                          mods=[ModCondition('phosphorylation', 'T', '183'),
                                ModCondition('phosphorylation', 'Y', '185')],
                          db_refs={'UP': 'P28482'})
    mapk3_invalid = Agent('MAPK3',
                          mods=[ModCondition('phosphorylation', 'T', '201'),
                                ModCondition('phosphorylation', 'Y', '203')],
                          db_refs={'UP': 'P27361'})

    st1 = Complex([mapk1_invalid, mapk3_invalid])
    res = sm.map_sites([st1])
    assert len(res) == 2
    valid_stmts = res[0]
    mapped_stmts = res[1]
    assert isinstance(valid_stmts, list)
    assert isinstance(mapped_stmts, list)
    assert len(valid_stmts) == 0
    assert len(mapped_stmts) == 1
    mapped_stmt = mapped_stmts[0]
    assert isinstance(mapped_stmt, MappedStatement)
    assert mapped_stmt.original_stmt == st1
    assert isinstance(mapped_stmt.mapped_mods, list)
    assert len(mapped_stmt.mapped_mods) == 4
    ms = mapped_stmt.mapped_stmt
    assert isinstance(ms, Statement)
    members = ms.members
    assert len(members) == 2
    agent1 = members[0]
    agent2 = members[1]
    assert agent1.name == 'MAPK1'
    assert len(agent1.mods) == 2
    assert agent1.mods[0].matches(ModCondition('phosphorylation', 'T', '185'))
    assert agent1.mods[1].matches(ModCondition('phosphorylation', 'Y', '187'))
    assert agent2.mods[0].matches(ModCondition('phosphorylation', 'T', '202'))
    assert agent2.mods[1].matches(ModCondition('phosphorylation', 'Y', '204'))

def test_site_map_modification():
    mapk1_invalid = Agent('MAPK1',
                          mods=[ModCondition('phosphorylation', 'T', '183'),
                                ModCondition('phosphorylation', 'Y', '185')],
                          db_refs={'UP': 'P28482'})
    mapk3_invalid = Agent('MAPK3',
                          mods=[ModCondition('phosphorylation', 'T', '201')],
                          db_refs={'UP': 'P27361'})
    map2k1_invalid = Agent('MAP2K1',
                           mods=[ModCondition('phosphorylation', 'S', '217'),
                                 ModCondition('phosphorylation', 'S', '221')],
                           db_refs={'UP': 'Q02750'})

    st1 = Phosphorylation(mapk1_invalid, mapk3_invalid, 'T', '203')
    st2 = Phosphorylation(map2k1_invalid, mapk1_invalid, 'Y', '218')

    res = sm.map_sites([st1, st2])

    assert len(res) == 2
    valid_stmts = res[0]
    mapped_stmts = res[1]
    assert isinstance(valid_stmts, list)
    assert isinstance(mapped_stmts, list)
    assert len(valid_stmts) == 0
    assert len(mapped_stmts) == 2
    # MAPK1 -> MAPK3
    mapped_stmt1 = mapped_stmts[0]
    assert isinstance(mapped_stmt1, MappedStatement)
    assert mapped_stmt1.original_stmt == st1
    assert isinstance(mapped_stmt1.mapped_mods, list)
    assert len(mapped_stmt1.mapped_mods) == 4 # FIXME
    ms = mapped_stmt1.mapped_stmt
    assert isinstance(ms, Statement)
    agent1 = ms.enz
    agent2 = ms.sub
    assert agent1.name == 'MAPK1'
    assert len(agent1.mods) == 2
    assert agent1.mods[0].matches(ModCondition('phosphorylation', 'T', '185'))
    assert agent1.mods[1].matches(ModCondition('phosphorylation', 'Y', '187'))
    assert agent2.mods[0].matches(ModCondition('phosphorylation', 'T', '202'))
    assert ms.residue == 'Y'
    assert ms.position == '204'

    # MAP2K1 -> MAPK1
    mapped_stmt2 = mapped_stmts[1]
    assert isinstance(mapped_stmt2, MappedStatement)
    assert mapped_stmt2.original_stmt == st2
    assert isinstance(mapped_stmt2.mapped_mods, list)
    assert len(mapped_stmt2.mapped_mods) == 5 # FIXME
    ms = mapped_stmt2.mapped_stmt
    assert isinstance(ms, Statement)
    agent1 = ms.enz
    agent2 = ms.sub
    assert agent1.name == 'MAP2K1'
    assert len(agent1.mods) == 2
    assert agent1.mods[0].matches(ModCondition('phosphorylation', 'S', '218'))
    assert agent1.mods[0].matches(ModCondition('phosphorylation', 'S', '222'))
    assert len(agent2.mods) == 2
    assert agent2.mods[0].matches(ModCondition('phosphorylation', 'T', '185'))
    assert agent2.mods[1].matches(ModCondition('phosphorylation', 'Y', '187'))
    # The incorrect phosphorylation residue is passed through to the new
    # statement unchanged
    assert ms.residue == 'Y'
    assert ms.position == '218'


"""
def test_site_map_activity_modification():
    mapk1_invalid = Agent('MAPK1', db_refs={'UP': 'P28482'})

    st1 = ActivityModification(mapk1_invalid,
                               ['PhosphorylationThreonine',
                                     'PhosphorylationTyrosine'],
                               ['183', '185'], 'increases',
                               'KinaseActivity')
    (valid, mapped) = sm.map_sites([st1])
    import ipdb; ipdb.set_trace()
"""

