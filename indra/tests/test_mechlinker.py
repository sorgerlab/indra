from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import networkx
from indra.statements import *
from indra.mechlinker import MechLinker, AgentState
from indra.mechlinker import _get_graph_reductions

def test_agent_state():
    mc = ModCondition('phosphorylation')
    mut = MutCondition('600', 'V', 'E')
    location = 'nucleus'
    bc = BoundCondition(Agent('x'), True)
    a = Agent('a', mods=[mc], mutations=[mut], bound_conditions=[bc],
              location=location)
    agent_state = AgentState(a)
    assert(agent_state.mods)
    assert(agent_state.mutations)
    assert(agent_state.bound_conditions)
    assert(agent_state.location)


def test_act_phos_to_af():
    act_st = Activation(Agent('A', activity=ActivityCondition('kinase', True)),
                        Agent('B'))
    phos_st = Phosphorylation(Agent('A'), Agent('B'))
    ml = MechLinker([act_st, phos_st])
    linked_stmts = ml.infer_active_forms(ml.statements)
    assert(len(linked_stmts) == 1)

def test_act_af_to_phos():
    act_st = Activation(Agent('A', activity=ActivityCondition('kinase', True)),
                        Agent('B'))
    af_st = ActiveForm(Agent('B', mods=[ModCondition('phosphorylation',
                                                     None, None, True)]),
                        'activity', True)
    ml = MechLinker([act_st, af_st])
    linked_stmts = ml.infer_modifications(ml.statements)
    assert(len(linked_stmts) == 1)

def test_reduce_activity_types():
    a1 = Agent('a', location='cytoplasm')
    a2 = Agent('a', location='nucleus')
    af1 = ActiveForm(a1, 'activity', True)
    af2 = ActiveForm(a2, 'kinase', True)
    af3 = ActiveForm(a1, 'catalytic',True)
    ml = MechLinker([af1, af2, af3])
    ml.gather_explicit_activities()
    ml.reduce_activities()
    assert(af1.activity == 'kinase')
    assert(af2.activity == 'kinase')
    assert(af3.activity == 'kinase')

def test_graph_reductions():
    G = networkx.DiGraph([('activity', 'kinase'),
                          ('catalytic', 'kinase'),
                          ('activity', 'catalytic'),
                          ('activity', 'phosphatase'),
                          ('catalytic', 'phosphatase')])
    reductions = _get_graph_reductions(G)
    assert(reductions == {'activity': 'catalytic',
                          'kinase': 'kinase',
                          'phosphatase': 'phosphatase',
                          'catalytic': 'catalytic'})
    G = networkx.DiGraph([('activity', 'kinase'),
                          ('catalytic', 'kinase'),
                          ('activity', 'catalytic')])
    reductions = _get_graph_reductions(G)
    assert(reductions == {'activity': 'kinase',
                          'catalytic': 'kinase',
                          'kinase': 'kinase'})
    G = networkx.DiGraph([('activity', 'kinase'),
                          ('catalytic', 'kinase'),
                          ('activity', 'catalytic'),
                          ('activity', 'transcription')])
    reductions = _get_graph_reductions(G)
    assert(reductions == {'activity': 'activity',
                          'transcription': 'transcription',
                          'catalytic': 'kinase',
                          'kinase': 'kinase'})

def test_base_agent():
    af = ActiveForm(Agent('a', mods=[ModCondition('phosphorylation')]),
                    'activity', True)
    ml = MechLinker([af])
    ml.gather_explicit_activities()

def test_require_active_forms_mod1():
    af = ActiveForm(Agent('a', mods=[ModCondition('phosphorylation')]),
                    'activity', True)
    ph = Phosphorylation(Agent('a'), Agent('b'))
    ml = MechLinker([af, ph])
    ml.gather_explicit_activities()
    ml.require_active_forms()
    assert(len(ml.statements) == 2)
    assert(ml.statements[1].enz.mods)

def test_require_active_forms_mod2():
    af = ActiveForm(Agent('a', mods=[ModCondition('phosphorylation')]),
                    'activity', True)
    af2 = ActiveForm(Agent('a', location='nucleus'), 'activity', True)
    ph = Phosphorylation(Agent('a'), Agent('b'))
    ml = MechLinker([af, af2, ph])
    ml.gather_explicit_activities()
    ml.require_active_forms()
    assert(len(ml.statements) == 4)
    assert(ml.statements[3].enz.location)

def test_require_active_forms_mod3():
    mc1 = ModCondition('phosphorylation', 'T', '185')
    mc2 = ModCondition('phosphorylation', 'Y', '187')
    af = ActiveForm(Agent('a', mods=[mc1, mc2]),
                    'kinase', True)
    ph = Phosphorylation(Agent('a'), Agent('b'))
    ml = MechLinker([af, ph])
    ml.gather_explicit_activities()
    ml.require_active_forms()
    assert(len(ml.statements) == 2)
    assert(len(ml.statements[1].enz.mods) == 2)

def test_require_active_forms_mod4():
    mc1 = ModCondition('phosphorylation', 'T', '185')
    mc2 = ModCondition('phosphorylation', 'Y', '187')
    af = ActiveForm(Agent('a', mods=[mc1, mc2]),
                    'kinase', True)
    ph = Phosphorylation(Agent('a', mods=[mc1]), Agent('b'))
    ml = MechLinker([af, ph])
    ml.gather_explicit_activities()
    ml.require_active_forms()
    assert(len(ml.statements) == 2)
    assert(len(ml.statements[1].enz.mods) == 2)

def test_require_active_forms_mod5():
    mc1 = ModCondition('phosphorylation', 'T', '185')
    mc2 = ModCondition('phosphorylation', 'Y', '187')
    mc3 = ModCondition('phosphorylation', 'S', '999')
    af = ActiveForm(Agent('a', mods=[mc1, mc2]),
                    'kinase', True)
    af2 = ActiveForm(Agent('a', mods=[mc3]),
                     'kinase', False)
    ph = Phosphorylation(Agent('a'), Agent('b'))
    ml = MechLinker([af, af2, ph])
    ml.gather_explicit_activities()
    ml.require_active_forms()
    assert(len(ml.statements) == 3)
    assert(len(ml.statements[2].enz.mods) == 2)

def test_require_active_forms_act1():
    af = ActiveForm(Agent('a', mods=[ModCondition('phosphorylation')]),
                    'activity', True)
    act = Activation(Agent('a'), Agent('b'))
    ml = MechLinker([af, act])
    ml.gather_explicit_activities()
    ml.require_active_forms()
    assert(len(ml.statements) == 2)
    assert(ml.statements[1].subj.mods)

def test_infer_activations():
    af = ActiveForm(Agent('a', mods=[ModCondition('phosphorylation')]),
                    'activity', True)
    phos = Phosphorylation(Agent('b'), Agent('a'))
    linked_stmts = MechLinker.infer_activations([af, phos])
    assert(len(linked_stmts) == 1)
    print(linked_stmts)

def test_replace_activations():
    af = ActiveForm(Agent('a', mods=[ModCondition('phosphorylation')]),
                    'activity', True)
    phos = Phosphorylation(Agent('b'), Agent('a'))
    act = Activation(Agent('b'), Agent('a'))
    ml = MechLinker([af, phos, act])
    ml.replace_activations()
    assert(len(ml.statements) == 2)
    print(ml.statements)

def test_infer_complexes():
    phos = Phosphorylation(Agent('b'), Agent('a'))
    linked_stmts = MechLinker.infer_complexes([phos])
    assert(len(linked_stmts) == 1)
    print(linked_stmts)

def test_replace_complexes():
    phos = Phosphorylation(Agent('b'), Agent('a'))
    cplx = Complex([Agent('a'), Agent('b')])
    ml = MechLinker([phos, cplx])
    ml.replace_complexes()
    assert(len(ml.statements) == 1)
    print(ml.statements)

def test_reduce_mods1():
    phos1 = Phosphorylation(Agent('b'), Agent('a'))
    phos2 = Phosphorylation(Agent('c'), Agent('a'), 'T')
    phos3 = Phosphorylation(Agent('d'), Agent('a'), 'T', '143')
    ml = MechLinker([phos1, phos2, phos3])
    ml.gather_modifications()
    ml.reduce_modifications()
    assert(len(ml.statements) == 3)
    for st in ml.statements:
        assert(st.residue == 'T')
        assert(st.position == '143')

def test_reduce_mods2():
    mc1 = ModCondition('phosphorylation', 'S', '123', False)
    mc2 = ModCondition('phosphorylation', 'S', None, True)
    mc3 = ModCondition('phosphorylation', 'T')
    mc4 = ModCondition('phosphorylation', 'T', '111')
    mc5 = ModCondition('phosphorylation', 'T', '999')
    mc6 = ModCondition('phosphorylation')
    mc7 = ModCondition('phosphorylation', None, '999')
    st1 = Activation(Agent('KRAS'), Agent('BRAF', mods=[mc1]))
    st2 = Activation(Agent('KRAS'), Agent('BRAF', mods=[mc2]))
    st3 = Activation(Agent('KRAS'), Agent('BRAF', mods=[mc3]))
    st4 = Activation(Agent('KRAS'), Agent('BRAF', mods=[mc4]))
    st5 = Activation(Agent('KRAS'), Agent('BRAF', mods=[mc5]))
    st6 = Activation(Agent('KRAS'), Agent('BRAF', mods=[mc6]))
    st7 = Activation(Agent('KRAS'), Agent('BRAF', mods=[mc7]))
    ml = MechLinker([st1, st2, st3, st4, st5, st6, st7])
    ml.gather_modifications()
    ml.reduce_modifications()
    assert(len(ml.statements) == 7)
    mc_red1 = ml.statements[0].obj.mods[0]
    mc_red2 = ml.statements[1].obj.mods[0]
    mc_red3 = ml.statements[2].obj.mods[0]
    mc_red4 = ml.statements[3].obj.mods[0]
    mc_red5 = ml.statements[4].obj.mods[0]
    mc_red6 = ml.statements[5].obj.mods[0]
    mc_red7 = ml.statements[6].obj.mods[0]
    # These ones stay the same because they shouldn't be reduced
    assert(mc_red1.__dict__ == mc1.__dict__)
    assert(mc_red3.__dict__ == mc3.__dict__)
    assert(mc_red4.__dict__ == mc4.__dict__)
    assert(mc_red5.__dict__ == mc5.__dict__)
    assert(mc_red6.__dict__ == mc6.__dict__)
    # mc2 has to be reduced to have position '123'
    assert(mc_red2.mod_type == 'phosphorylation')
    assert(mc_red2.residue == 'S')
    assert(mc_red2.position == '123')
    assert(mc_red2.is_modified == True)
    # mc7 has to be reduced to have residue 'T'
    assert(mc_red7.mod_type == 'phosphorylation')
    assert(mc_red7.residue == 'T')
    assert(mc_red7.position == '999')
    assert(mc_red7.is_modified == True)
