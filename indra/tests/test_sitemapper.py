from protmapper import MappedSite
from indra.statements import *
from indra.util import unicode_strs
from indra.preassembler.sitemapper import default_mapper as sm, MappedStatement
from indra.preassembler.sitemapper import _valid_position_str


def test_check_agent_mod():
    mapk1_valid = Agent('MAPK1',
                        mods=[ModCondition('phosphorylation', 'T', '185'),
                              ModCondition('phosphorylation', 'Y', '187')],
                        db_refs={'UP': 'P28482'})
    mapped_sites_valid, _ = sm._map_agent_sites(mapk1_valid)
    assert not mapped_sites_valid, mapped_sites_valid

    mapk1_invalid = Agent('MAPK1',
                          mods=[ModCondition('phosphorylation', 'T', '183'),
                                ModCondition('phosphorylation', 'Y', '185')],
                          db_refs={'UP': 'P28482'})
    mapped_sites_invalid, new_agent = sm._map_agent_sites(mapk1_invalid)
    assert len(mapped_sites_invalid) == 2
    assert isinstance(mapped_sites_invalid[0], MappedSite)
    map183 = mapped_sites_invalid[0]
    assert (map183.up_id, map183.orig_res, map183.orig_pos,
            map183.mapped_res, map183.mapped_pos) == \
        ('P28482', 'T', '183', 'T', '185'), map183
    map185 = mapped_sites_invalid[1]
    assert (map185.up_id, map185.orig_res, map185.orig_pos,
            map185.mapped_res, map185.mapped_pos) == \
        ('P28482', 'Y', '185', 'Y', '187'), map183
    assert len(new_agent.mods) == 2
    assert new_agent.mods[0].matches(ModCondition('phosphorylation',
                                                  'T', '185'))
    assert new_agent.mods[1].matches(ModCondition('phosphorylation',
                                                  'Y', '187'))
    assert unicode_strs((mapk1_valid, mapk1_invalid, mapped_sites_invalid,
                         mapped_sites_valid, map183, map185, new_agent))


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

    st1 = Phosphorylation(mapk1_invalid, mapk3_invalid, 'Y', '203')
    st2 = Phosphorylation(map2k1_invalid, mapk1_invalid, 'Y', '218')
    res = sm.map_sites([st1, st2])

    assert len(res) == 2
    valid_stmts, mapped_stmts = res
    assert isinstance(valid_stmts, list)
    assert isinstance(mapped_stmts, list)
    assert len(valid_stmts) == 0
    assert len(mapped_stmts) == 2
    # MAPK1 -> MAPK3
    mapped_stmt1 = mapped_stmts[0]
    assert isinstance(mapped_stmt1, MappedStatement)
    assert mapped_stmt1.original_stmt == st1
    assert isinstance(mapped_stmt1.mapped_mods, list)
    assert len(mapped_stmt1.mapped_mods) == 4, mapped_stmt1.mapped_mods
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
    assert len(mapped_stmt2.mapped_mods) == 5, mapped_stmt2.mapped_mods
    ms = mapped_stmt2.mapped_stmt
    assert isinstance(ms, Statement)
    agent1 = ms.enz
    agent2 = ms.sub
    assert agent1.name == 'MAP2K1'
    assert len(agent1.mods) == 2
    assert agent1.mods[0].matches(ModCondition('phosphorylation', 'S', '218'))
    assert agent1.mods[1].matches(ModCondition('phosphorylation', 'S', '222'))
    assert len(agent2.mods) == 2
    assert agent2.mods[0].matches(ModCondition('phosphorylation', 'T', '185'))
    assert agent2.mods[1].matches(ModCondition('phosphorylation', 'Y', '187'))
    # The incorrect phosphorylation residue is passed through to the new
    # statement unchanged
    assert ms.residue == 'Y'
    assert ms.position == '218'
    # Check for unicode
    assert unicode_strs((mapk1_invalid, mapk3_invalid, map2k1_invalid, st1,
                         st2, res, valid_stmts, mapped_stmts))


def test_site_map_activity_modification():
    mc = [ModCondition('phosphorylation', 'T', '183'),
          ModCondition('phosphorylation', 'Y', '185')]
    mapk1 = Agent('MAPK1', mods=mc, db_refs={'UP': 'P28482'})

    st1 = ActiveForm(mapk1, 'kinase', True)
    (valid, mapped) = sm.map_sites([st1])
    assert len(valid) == 0
    assert len(mapped) == 1
    ms = mapped[0]
    mm = ms.mapped_mods
    assert (mm[0].gene_name, mm[0].orig_res, mm[0].orig_pos, mm[0].mapped_res,
            mm[0].mapped_pos) == ('MAPK1', 'T', '183', 'T', '185')
    assert (mm[1].gene_name, mm[1].orig_res, mm[1].orig_pos, mm[1].mapped_res,
            mm[1].mapped_pos) == ('MAPK1', 'Y', '185', 'Y', '187')
    assert ms.original_stmt == st1
    assert ms.mapped_stmt.agent.mods[0].matches(ModCondition('phosphorylation',
                                                             'T', '185'))
    assert ms.mapped_stmt.agent.mods[1].matches(ModCondition('phosphorylation',
                                                             'Y', '187'))
    assert unicode_strs((mc, mapk1, st1, valid, mapped))


def test_site_map_selfmodification():
    mapk1_invalid = Agent('MAPK1',
                          mods=[ModCondition('phosphorylation', 'T', '183')],
                          db_refs={'UP': 'P28482'})
    st1 = Autophosphorylation(mapk1_invalid, 'Y', '185')
    (valid, mapped) = sm.map_sites([st1])
    assert len(valid) == 0
    assert len(mapped) == 1
    mapped_stmt = mapped[0]
    mm = mapped_stmt.mapped_mods
    assert (mm[0].gene_name, mm[0].orig_res, mm[0].orig_pos, mm[0].mapped_res,
            mm[0].mapped_pos) == ('MAPK1', 'T', '183', 'T', '185')
    assert (mm[1].gene_name, mm[1].orig_res, mm[1].orig_pos, mm[1].mapped_res,
            mm[1].mapped_pos) == ('MAPK1', 'Y', '185', 'Y', '187')
    ms = mapped_stmt.mapped_stmt
    agent1 = ms.enz
    assert agent1.mods[0].matches(ModCondition('phosphorylation', 'T', '185'))
    assert ms.residue == 'Y'
    assert ms.position == '187'
    assert unicode_strs((mapk1_invalid, st1, valid, mapped))


# The following Statements are all handled by the same block of code and hence
# can be tested in similar fashion


def test_site_map_complex():
    (mapk1_invalid, mapk3_invalid) = get_invalid_mapks()
    st1 = Complex([mapk1_invalid, mapk3_invalid])
    res = sm.map_sites([st1])
    check_validated_mapks(res, st1)


def test_site_map_gef():
    (mapk1_invalid, mapk3_invalid) = get_invalid_mapks()
    st1 = Gef(mapk1_invalid, mapk3_invalid)
    res = sm.map_sites([st1])
    check_validated_mapks(res, st1)


def test_site_map_gap():
    (mapk1_invalid, mapk3_invalid) = get_invalid_mapks()
    st1 = Gap(mapk1_invalid, mapk3_invalid)
    res = sm.map_sites([st1])
    check_validated_mapks(res, st1)


def test_site_map_activation():
    (mapk1_invalid, mapk3_invalid) = get_invalid_mapks()
    st1 = Activation(mapk1_invalid, mapk3_invalid, 'kinase')
    res = sm.map_sites([st1])
    check_validated_mapks(res, st1)


def test_site_map_hgnc():
    """Make sure site mapping is done even if only HGNC ID is given."""
    (mapk1_invalid, mapk3_invalid) = get_invalid_mapks()
    mapk1_invalid.db_refs = {'HGNC': '6871'}
    st1 = ActiveForm(mapk1_invalid, 'kinase', True)
    (valid, mapped) = sm.map_sites([st1])
    assert len(valid) == 0
    assert len(mapped) == 1


def test_site_map_within_bound_condition():
    # Here, we test to make sure that agents within a bound condition are
    # site-mapped
    (mapk1_invalid, mapk3_invalid) = get_invalid_mapks()

    # Add an agent to the bound condition for the object of the statement
    mapk3_invalid.bound_conditions = [BoundCondition(mapk1_invalid)]
    st1 = Activation(mapk1_invalid, mapk3_invalid, 'kinase')

    # Map sites
    res = sm.map_sites([st1])

    # Extract the mapped statement
    mapped_statements = res[1]
    assert len(mapped_statements) == 1
    mapped_s = mapped_statements[0].mapped_stmt

    # Verify that the agent in the object's bound condition got site-mapped
    validate_mapk1(mapped_s.obj.bound_conditions[0].agent)


def test_invalid_position():
    stmt = Phosphorylation._from_json({
        'enz': {'name': 'CFD'},
        'sub': {'name': 'HP'},
        'residue': 'F',
        'position': '2.59'
    })
    valid, mapped = sm.map_sites(stmts=[stmt])
    assert not valid
    assert not mapped


def test_valid_pos_str():
    assert _valid_position_str('5') is True
    assert _valid_position_str('005') is False
    assert _valid_position_str('5.5') is False
    assert _valid_position_str('xxx') is False
    assert _valid_position_str(None) is True
    assert _valid_position_str('') is False


def get_invalid_mapks():
    """A handy function for getting the invalid MAPK agents we want."""
    mapk1_invalid = Agent('MAPK1',
                          mods=[ModCondition('phosphorylation', 'T', '183'),
                                ModCondition('phosphorylation', 'Y', '185')],
                          db_refs={'UP': 'P28482'})
    mapk3_invalid = Agent('MAPK3',
                          mods=[ModCondition('phosphorylation', 'T', '201'),
                                ModCondition('phosphorylation', 'Y', '203')],
                          db_refs={'UP': 'P27361'})
    assert unicode_strs((mapk1_invalid, mapk3_invalid))
    return (mapk1_invalid, mapk3_invalid)


def check_validated_mapks(res, st1):
    """Validate that the invalid MAPKs have been fixed appropriately."""
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
    agents = ms.agent_list()
    assert len(agents) == 2
    agent1 = agents[0]
    agent2 = agents[1]
    validate_mapk1(agent1)
    assert agent2.mods[0].matches(ModCondition('phosphorylation', 'T', '202'))
    assert agent2.mods[1].matches(ModCondition('phosphorylation', 'Y', '204'))
    assert unicode_strs((res, st1))


def validate_mapk1(agent1):
    assert agent1.name == 'MAPK1'
    assert len(agent1.mods) == 2
    assert agent1.mods[0].matches(ModCondition('phosphorylation', 'T', '185'))
    assert agent1.mods[1].matches(ModCondition('phosphorylation', 'Y', '187'))

