from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.preassembler.grounding_mapper import GroundingMapper, \
                                                default_grounding_map

def get_statements():
    statements = []

    egf = Agent('EGF')
    egfr = Agent('EGFR')
    st = Complex([egf, egfr])
    statements.append(st)

    egfre = Agent('EGFR', bound_conditions=[BoundCondition(egf, True)])
    egfre = Agent('EGFR', bound_conditions=[BoundCondition(egf, True)])
    st = Complex([egfre, egfre])
    statements.append(st)

    egfrdimer = Agent('EGFR', bound_conditions=[BoundCondition(egfr, True)])
    st = Transphosphorylation(egfrdimer, 'Y')
    statements.append(st)

    egfrpY = Agent('EGFR', mods=[ModCondition('phosphorylation', 'Y')])
    grb2 = Agent('GRB2')
    st = Complex([egfrpY, grb2])
    statements.append(st)


    grb2bound = Agent('GRB2', bound_conditions=[BoundCondition(egfr, True)])
    sos1 = Agent('SOS1')
    st = Complex([grb2bound, sos1])
    statements.append(st)

    hras = Agent('HRAS')
    kras = Agent('KRAS')
    nras = Agent('NRAS')
    gdp = Agent('GDP')
    for ras in [hras, kras, nras]:
        st = Complex([ras, gdp])
        statements.append(st)

    sos1bound = Agent('SOS1', bound_conditions=[BoundCondition(grb2, True)])
    hras_gdp = Agent('HRAS', bound_conditions=[BoundCondition(gdp, True)])
    kras_gdp = Agent('KRAS', bound_conditions=[BoundCondition(gdp, True)])
    nras_gdp = Agent('NRAS', bound_conditions=[BoundCondition(gdp, True)])
    for ras_gdp in [hras_gdp, kras_gdp, nras_gdp]:
        st = Complex([sos1bound, ras_gdp])
        statements.append(st)
        st = ActiveForm(ras_gdp, 'activity', False)
        statements.append(st)


    hras_bound = Agent('HRAS', bound_conditions=[BoundCondition(sos1, True)])
    kras_bound = Agent('KRAS', bound_conditions=[BoundCondition(sos1, True)])
    nras_bound = Agent('NRAS', bound_conditions=[BoundCondition(sos1, True)])
    sos1bound = Agent('SOS1', bound_conditions=[BoundCondition(grb2, True)])
    for ras_bound in [hras_bound, kras_bound, nras_bound]:
        st = Complex([sos1bound, ras_bound])
        statements.append(st)


    gtp = Agent('GTP')
    hras_gtp = Agent('HRAS', bound_conditions=[BoundCondition(gtp, True)])
    kras_gtp = Agent('KRAS', bound_conditions=[BoundCondition(gtp, True)])
    nras_gtp = Agent('NRAS', bound_conditions=[BoundCondition(gtp, True)])
    braf = Agent('BRAF')
    for ras_gtp in [hras_gtp, kras_gtp, nras_gtp]:
        st = Complex([ras_gtp, braf])
        statements.append(st)
        st = ActiveForm(ras_gtp, 'activity', True)
        statements.append(st)

    hras_braf = Agent('BRAF', bound_conditions=[BoundCondition(hras, True)])
    kras_braf = Agent('BRAF', bound_conditions=[BoundCondition(kras, True)])
    nras_braf = Agent('BRAF', bound_conditions=[BoundCondition(nras, True)])
    for braf1 in [hras_braf, kras_braf, nras_braf]:
        for braf2 in [hras_braf, kras_braf, nras_braf]:
            st = Complex([braf1, braf2])
            statements.append(st)

    braf_bound = Agent('BRAF', bound_conditions=[BoundCondition(braf, True)])
    st = Transphosphorylation(braf_bound)
    statements.append(st)

    braf_phos = Agent('BRAF', mods=[ModCondition('phosphorylation')])
    mek1 = Agent('MAP2K1')
    mek2 = Agent('MAP2K2')
    st = ActiveForm(braf_phos, 'kinase', True)
    statements.append(st)
    st = Phosphorylation(braf_phos, mek1)
    statements.append(st)
    st = Phosphorylation(braf_phos, mek2)
    statements.append(st)

    mek1_phos = Agent('MAP2K1', mods=[ModCondition('phosphorylation')])
    mek2_phos = Agent('MAP2K2', mods=[ModCondition('phosphorylation')])
    mapk1 = Agent('MAPK1')
    mapk3 = Agent('MAPK3')
    st = ActiveForm(mek1_phos, 'kinase', True)
    statements.append(st)
    st = ActiveForm(mek2_phos, 'kinase', True)
    statements.append(st)
    st = Phosphorylation(braf_phos, mek1)
    statements.append(st)
    st = Phosphorylation(braf_phos, mek2)
    statements.append(st)
    for mek in [mek1_phos, mek2_phos]:
        for erk in [mapk1, mapk3]:
            st = Phosphorylation(mek, erk)

    for st in statements:
        st.belief = 1
        st.evidence.append(Evidence(source_api='assertion'))

    # Update the statements with grounding info. To do this, we set the "text"
    # field of the db_refs to copy from the agent name, then run the grounding
    # mapper
    for st in statements:
        for ag in st.agent_list():
            if ag is None:
                continue
            else:
                ag.db_refs = {'TEXT': ag.name}
    # Now load the grounding map and run
    gm = GroundingMapper(default_grounding_map)
    mapped_stmts = gm.map_agents(statements)
    # This shouldn't change anything, but just in case...
    renamed_stmts = gm.rename_agents(mapped_stmts)
    return renamed_stmts


