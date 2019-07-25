import pandas as pd
from indra.assemblers.indranet import IndraNetAssembler
from indra.statements import *


ev1 = Evidence(pmid='1')
ev2 = Evidence(pmid='2')
ev3 = Evidence(pmid='3')
st1 = Activation(Agent('a', db_refs={'HGNC': 1}), Agent('b'), evidence=[ev1])
st2 = Inhibition(Agent('a', db_refs={'HGNC': 1}), Agent('c'),
                 evidence=[ev1, ev2, ev3])
st2.belief = 0.76
st3 = Activation(Agent('b'), Agent('d'))
st4 = ActiveForm(Agent('e'), None, True)  # 1 agent
st5 = Complex([Agent('c'), Agent('f'), Agent('g')])
st6 = Complex([Agent('h'), Agent('i'), Agent('j'), Agent('b')])
st7 = Phosphorylation(None, Agent('x'))


def test_simple_assembly():
    ia = IndraNetAssembler([st1, st2, st3, st4, st5, st6, st7])
    g = ia.make_model()
    assert len(g.nodes) == 6
    assert len(g.edges) == 9
    print(g.edges)
    # Stmt with 1 agent should not be added
    assert 'e' not in g.nodes
    # Complex with more than 3 agents should not be added
    assert ('f', 'g', 0) in g.edges
    assert ('h', 'i', 0) not in g.edges
    # Test node attributes
    assert g.nodes['a']['ns'] == 'HGNC', g.nodes['a']['ns']
    assert g.nodes['a']['id'] == 1
    # Test edge attributes
    e = g['a']['c'][0]
    assert e['stmt_type'] == 'Inhibition'
    assert e['belief'] == 0.76
    assert e['evidence_count'] == 3
    assert g['b']['d'][0]['evidence_count'] == 0


def test_exclude_stmts():
    ia = IndraNetAssembler([st1, st2, st3])
    g = ia.make_model(exclude_stmts=['Inhibition'])
    assert len(g.nodes) == 3
    assert len(g.edges) == 2
    assert 'c' not in g.nodes
    assert ('a', 'c', 0) not in g.edges


def test_complex_members():
    ia = IndraNetAssembler([st1, st6])
    g = ia.make_model(complex_members=4)
    assert len(g.nodes) == 5
    assert len(g.edges) == 13, len(g.edges)
    assert ('h', 'i', 0) in g.edges
    assert ('i', 'h', 0) in g.edges


def test_make_df():
    ia = IndraNetAssembler([st1, st2, st3, st4, st5, st6])
    df = ia.make_df()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 9
    assert set(df.columns) == {
        'agA_name', 'agB_name', 'agA_ns', 'agA_id', 'agB_ns', 'agB_id',
        'stmt_type', 'evidence_count', 'stmt_hash', 'belief'}
