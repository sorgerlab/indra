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

