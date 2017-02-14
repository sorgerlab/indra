from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import networkx
from indra.statements import *
from indra.mechlinker import MechLinker
from indra.mechlinker import get_graph_reductions

def test_act_phos_to_af():
    act_st = Activation(Agent('A', activity=ActivityCondition('kinase', True)),
                        Agent('B'))
    phos_st = Phosphorylation(Agent('A'), Agent('B'))
    ml = MechLinker([act_st, phos_st])
    linked_stmts = ml.link_statements()
    assert(len(linked_stmts) == 1)

def test_act_af_to_phos():
    act_st = Activation(Agent('A', activity=ActivityCondition('kinase', True)),
                        Agent('B'))
    af_st = ActiveForm(Agent('B', mods=[ModCondition('phosphorylation',
                                                     None, None, True)]),
                        'activity', True)
    ml = MechLinker([act_st, af_st])
    linked_stmts = ml.link_statements()
    assert(len(linked_stmts) == 1)

def test_reduce_activity_types():
    a1 = Agent('a', location='cytoplasm')
    a2 = Agent('a', location='nucleus')
    af1 = ActiveForm(a1, 'activity', True)
    af2 = ActiveForm(a2, 'kinase', True)
    af3 = ActiveForm(a1, 'catalytic',True)
    ml = MechLinker([af1, af2, af3])
    ml.get_explicit_activities()
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
    reductions = get_graph_reductions(G)
    assert(reductions == {'activity': 'catalytic',
                          'kinase': 'kinase',
                          'phosphatase': 'phosphatase',
                          'catalytic': 'catalytic'})
    G = networkx.DiGraph([('activity', 'kinase'),
                          ('catalytic', 'kinase'),
                          ('activity', 'catalytic')])
    reductions = get_graph_reductions(G)
    assert(reductions == {'activity': 'kinase',
                          'catalytic': 'kinase',
                          'kinase': 'kinase'})
    G = networkx.DiGraph([('activity', 'kinase'),
                          ('catalytic', 'kinase'),
                          ('activity', 'catalytic'),
                          ('activity', 'transcription')])
    reductions = get_graph_reductions(G)
    assert(reductions == {'activity': 'activity',
                          'transcription': 'transcription',
                          'catalytic': 'kinase',
                          'kinase': 'kinase'})
