from indra.preassembler.sitemapper import default_mapper as sm
from indra.statements import *

def test_check_agent_mod():
    mapk1_valid = Agent('MAPK1',
                         mods=['PhosphorylationThreonine',
                                'PhosphorylationTyrosine'],
                          mod_sites=['185', '187'], db_refs={'UP': 'P28482'})
    res_valid = sm.check_agent_mod(mapk1_valid)
    assert len(res_valid) == 2
    assert res_valid[0] == {}
    assert isinstance(res_valid[1], Agent)
    assert res_valid[1].matches(mapk1_valid)

    mapk1_invalid = Agent('MAPK1',
                          mods=['PhosphorylationThreonine',
                                'PhosphorylationTyrosine'],
                          mod_sites=['183', '185'], db_refs={'UP': 'P28482'})
    res_invalid = sm.check_agent_mod(mapk1_invalid)
    assert len(res_invalid) == 2
    assert isinstance(res_invalid[0], dict)
    assert isinstance(res_invalid[1], Agent)
    invalid_sites = res_invalid[0]
    assert len(invalid_sites.keys()) == 2
    map183 = invalid_sites[('MAPK1', 'T', '183')]
    assert len(map183) == 3
    assert map183[0] == 'T'
    assert map183[1] == '185'
    map185 = invalid_sites[('MAPK1', 'Y', '185')]
    assert len(map185) == 3
    assert map185[0] == 'Y'
    assert map185[1] == '187'
    new_agent = res_invalid[1]
    assert new_agent.mods == ['PhosphorylationThreonine',
                             'PhosphorylationTyrosine']
    assert new_agent.mod_sites == ['185', '187']

"""
def test_site_map_complex():
    mapk1_invalid = Agent('MAPK1',
                          mods=['PhosphorylationThreonine',
                                'PhosphorylationThreonine'],
                          mod_sites=['183', '185'], db_refs={'UP': 'P28482'})
    mapk3_invalid = Agent('MAPK3',
                            mods=['PhosphorylationThreonine',
                                  'PhosphorylationThreonine'],
                            mod_sites=['201', '203'])

    st1 = Complex([mapk1_invalid, mapk3_invalid])
    res = sm.check_(st1)
    # Also check case where statement has site that is incorrect but doesn't map
"""
