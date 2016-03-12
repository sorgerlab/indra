from indra.preassembler.site_mapper import default_mapper as sm
from indra.statements import *

def test_check_agent_mod():
    mapk1_valid = Agent('MAPK1', mods=['PhosphorylationThreonine'],
                                 mod_sites=['185'],
                                 db_refs={'UP': 'P28482'})
    res_valid = sm.check_agent_mod(mapk1_valid)
    assert res_valid == []

    mapk1_invalid = Agent('MAPK1', mods=['PhosphorylationThreonine'],
                                   mod_sites=['183'],
                                   db_refs={'UP': 'P28482'})
    import ipdb; ipdb.set_trace()
    res_invalid = sm.check_agent_mod(mapk1_invalid)
    assert res_invalid == [('MAPK1', 'T', '183')]


