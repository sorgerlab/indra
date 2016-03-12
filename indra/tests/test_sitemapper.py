from indra.preassembler.sitemapper import default_mapper as sm
from indra.statements import *

def test_check_agent_mod():
    mapk1_valid = Agent('MAPK1', mods=['PhosphorylationThreonine'],
                                 mod_sites=['185'],
                                 db_refs={'UP': 'P28482'})
    res_valid = sm.check_agent_mod(mapk1_valid)
    assert res_valid == {}

    mapk1_invalid = Agent('MAPK1', mods=['PhosphorylationThreonine'],
                                   mod_sites=['183'],
                                   db_refs={'UP': 'P28482'})
    res_invalid = sm.check_agent_mod(mapk1_invalid)
    import ipdb; ipdb.set_trace()
    assert res_invalid.keys()[0] == ('MAPK1', 'T', '183')
    mapped_site = res_invalid.values()[0]
    assert mapped_site[0] == 'T'
    assert mapped_site[1] == '185'




