import numpy as np
import pandas as pd
import networkx as nx
from indra.statements import *
from indra.assemblers.indranet.net import default_sign_dict
from indra.assemblers.indranet import IndraNetAssembler, IndraNet


ev1 = Evidence(pmid='1')
ev2 = Evidence(pmid='2')
ev3 = Evidence(pmid='3')
st1 = Activation(Agent('a', db_refs={'HGNC': '1'}), Agent('b'), evidence=[ev1])
st2 = Inhibition(Agent('a', db_refs={'HGNC': '1'}), Agent('c'),
                 evidence=[ev1, ev2, ev3])
st2.belief = 0.76
st3 = Activation(Agent('b'), Agent('d'))
st4 = ActiveForm(Agent('e'), None, True)  # 1 agent
st5 = Complex([Agent('c'), Agent('f'), Agent('g')])
st6 = Complex([Agent('h'), Agent('i'), Agent('j'), Agent('b')])
st7 = Phosphorylation(None, Agent('x'))
st8 = Conversion(Agent('PI3K'), [Agent('PIP2')], [Agent('PIP3')])


# Test assembly from assembler side
def test_simple_assembly():
    ia = IndraNetAssembler([st1, st2, st3, st4, st5, st6, st7])
    g = ia.make_model()
    assert len(g.nodes) == 6
    assert len(g.edges) == 9
    # Stmt with 1 agent should not be added
    assert 'e' not in g.nodes
    # Complex with more than 3 agents should not be added
    assert ('f', 'g', 0) in g.edges
    assert ('h', 'i', 0) not in g.edges
    # Test node attributes
    assert g.nodes['a']['ns'] == 'HGNC', g.nodes['a']['ns']
    assert g.nodes['a']['id'] == '1'
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
        'stmt_type', 'evidence_count', 'stmt_hash', 'belief', 'source_counts',
        'initial_sign'}


# Test assembly from IndraNet directly
def test_from_df():
    ia = IndraNetAssembler([st1, st2, st3, st4, st5, st6, st7])
    df = ia.make_df()
    net = IndraNet.from_df(df)
    assert len(net.nodes) == 6
    assert len(net.edges) == 9
    # Stmt with 1 agent should not be added
    assert 'e' not in net.nodes
    # Complex with more than 3 agents should not be added
    assert ('f', 'g', 0) in net.edges
    assert ('h', 'i', 0) not in net.edges
    # Test node attributes
    assert net.nodes['a']['ns'] == 'HGNC', net.nodes['a']['ns']
    assert net.nodes['a']['id'] == '1'
    # Test edge attributes
    e = net['a']['c'][0]
    assert e['stmt_type'] == 'Inhibition'
    assert e['belief'] == 0.76
    assert e['evidence_count'] == 3
    assert net['b']['d'][0]['evidence_count'] == 0


ab1 = Activation(Agent('a'), Agent('b'), evidence=[
    Evidence(source_api='sparser')])
ab2 = Phosphorylation(Agent('a'), Agent('b'),evidence=[
    Evidence(source_api='sparser'), Evidence(source_api='reach')])
ab3 = Inhibition(Agent('a'), Agent('b'), evidence=[
    Evidence(source_api='sparser'), Evidence(source_api='reach')])
ab4 = IncreaseAmount(Agent('a'), Agent('b'), evidence=[
    Evidence(source_api='trips')])
bc1 = Activation(Agent('b'), Agent('c'), evidence=[
    Evidence(source_api='trips')])
bc2 = Inhibition(Agent('b'), Agent('c'), evidence=[
    Evidence(source_api='trips'), Evidence(source_api='reach')])
bc3 = IncreaseAmount(Agent('b'), Agent('c'), evidence=[
    Evidence(source_api='sparser'), Evidence(source_api='reach')])
bc4 = DecreaseAmount(Agent('b'), Agent('c'), evidence=[
    Evidence(source_api='reach'), Evidence(source_api='trips')])


def test_to_digraph():
    ia = IndraNetAssembler([ab1, ab2, ab3, ab4, bc1, bc2, bc3, bc4])
    df = ia.make_df()
    net = IndraNet.from_df(df)
    assert len(net.nodes) == 3
    assert len(net.edges) == 8
    digraph = net.to_digraph(weight_mapping=_weight_mapping)
    assert len(digraph.nodes) == 3
    assert len(digraph.edges) == 2
    assert set([
        stmt['stmt_type'] for stmt in digraph['a']['b']['statements']]) == {
            'Activation', 'Phosphorylation', 'Inhibition', 'IncreaseAmount'}
    assert all(digraph.edges[e].get('belief', False) for e in digraph.edges)
    assert all(isinstance(digraph.edges[e]['belief'],
                          (float, np.longfloat)) for e in digraph.edges)
    assert all(digraph.edges[e].get('weight', False) for e in digraph.edges)
    assert all(isinstance(digraph.edges[e]['weight'],
                          (float, np.longfloat)) for e in digraph.edges)
    digraph_from_df = IndraNet.digraph_from_df(df)
    assert nx.is_isomorphic(digraph, digraph_from_df)


def test_to_signed_graph():
    ia = IndraNetAssembler([ab1, ab2, ab3, ab4, bc1, bc2, bc3, bc4])
    df = ia.make_df()
    net = IndraNet.from_df(df)
    signed_graph = net.to_signed_graph(
        sign_dict=default_sign_dict,
        weight_mapping=_weight_mapping)
    assert len(signed_graph.nodes) == 3
    assert len(signed_graph.edges) == 4
    assert set([stmt['stmt_type'] for stmt in
                signed_graph['a']['b'][0]['statements']]) == {
                    'Activation', 'IncreaseAmount'}
    assert set([stmt['stmt_type'] for stmt in
                signed_graph['a']['b'][1]['statements']]) == {'Inhibition'}
    assert set([stmt['stmt_type'] for stmt in
                signed_graph['b']['c'][0]['statements']]) == {
                    'Activation', 'IncreaseAmount'}
    assert set([stmt['stmt_type'] for stmt in
                signed_graph['b']['c'][1]['statements']]) == {
                    'Inhibition', 'DecreaseAmount'}
    assert all(signed_graph.edges[e].get('belief', False) for e in
               signed_graph.edges)
    assert all(isinstance(signed_graph.edges[e]['belief'],
                          (float, np.longfloat)) for e in signed_graph.edges)
    assert all(signed_graph.edges[e].get('weight', False) for e in
               signed_graph.edges)
    assert all(isinstance(signed_graph.edges[e]['weight'],
                          (float, np.longfloat)) for e in signed_graph.edges)


def _weight_mapping(G):
    for edge in G.edges:
        G.edges[edge]['weight'] = 1 - G.edges[edge]['belief']
    return G


def test_initial_signs():
    a = Event(Concept('a'), QualitativeDelta(polarity=1))
    b = Event(Concept('b'), QualitativeDelta(polarity=1))
    c = Event(Concept('c'), QualitativeDelta(polarity=-1))
    d = Event(Concept('d'), QualitativeDelta(polarity=-1))
    st1 = Influence(a, b)
    st2 = Influence(b, c)
    st3 = Influence(c, d)
    st4 = Influence(b, d)
    ia = IndraNetAssembler([st1, st2, st3, st4])
    sg = ia.make_model(graph_type='signed')
    assert len(sg.nodes) == 4
    assert len(sg.edges) == 4
    assert ('a', 'b', 0) in sg.edges
    assert ('b', 'c', 0) not in sg.edges
    assert ('b', 'c', 1) in sg.edges
    assert ('c', 'd', 0) in sg.edges
    assert ('c', 'd', 1) not in sg.edges
    assert ('b', 'd', 0) not in sg.edges
    assert ('b', 'd', 1) in sg.edges


def test_conversion():
    ia = IndraNetAssembler([st8])
    ug = ia.make_model(graph_type='multi_graph')
    assert len(ug.nodes) == 3
    assert len(ug.edges) == 2, ug.edges
    sg = ia.make_model(graph_type='signed')
    assert len(sg.nodes) == 3
    assert len(sg.edges) == 2, sg.edges
    assert ('PI3K', 'PIP3', 0) in sg.edges, sg.edges
    assert ('PI3K', 'PIP2', 1) in sg.edges, sg.edges
