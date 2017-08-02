from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers import CyJSAssembler

mek = Agent('MAP2K1', db_refs={'HGNC': '6840'})
erk = Agent('MAPK1', db_refs={'UP': 'P28482'})
dusp = Agent('DUSP4')
st_phos = Phosphorylation(mek, erk)
st_phos_Y = Phosphorylation(mek, erk, residue='Y')
st_phos_T = Phosphorylation(mek, erk, residue='T')
st_dephos = Dephosphorylation(dusp, erk)
st_complex = Complex([mek, erk, dusp])
st_act = Activation(mek, erk)
st_gef = Gef(Agent('SOS1'), Agent('HRAS'))
st_gap = Gap(Agent('RASA1'), Agent('HRAS'))
st_incamount = IncreaseAmount(Agent('TP53'), Agent('MDM2'))
st_decamount = DecreaseAmount(Agent('MDM2'), Agent('TP53'))
st_act2 = Inhibition(dusp, erk)
st_cited = Phosphorylation(mek, erk, evidence=Evidence(pmid='12345',
                                              text='MEK phosphorylates ERK'))
st_cited2 = Phosphorylation(mek, erk, evidence=Evidence(pmid='api35',
                                              text='MEK phosphorylates ERK'))
st_selfmod = Autophosphorylation(Agent('AKT1'), 'S', '473')

def test_act():
    cja = CyJSAssembler()
    cja.add_statements([st_act, st_act2])
    cja.make_model()
    assert(len(cja._nodes) == 3)
    assert(len(cja._edges) == 2)
    polarities = [edge['data']['polarity'] for edge in cja._edges]
    assert(len(set(polarities))==2)
    assert('positive' in polarities)
    assert('negative' in polarities)
    db_refs = [node['data']['db_refs'] for node in cja._nodes]
    for node in cja._nodes:
        if node['data']['name'] == 'MAP2K1':
            assert(node['data']['db_refs'].get('HGNC'))
        if node['data']['name'] == 'MAPK1':
            assert(node['data']['db_refs'].get('UniProt'))
        if node['data']['name'] == 'DUSP4':
            assert(not node['data']['db_refs'])

def test_regamount():
    cja = CyJSAssembler()
    cja.add_statements([st_incamount, st_decamount])
    cja.make_model()
    assert(len(cja._nodes) == 2)
    assert(len(cja._edges) == 2)
    polarities = [edge['data']['polarity'] for edge in cja._edges]
    assert(len(set(polarities))==2)
    assert('positive' in polarities)
    assert('negative' in polarities)

def test_ras():
    cja = CyJSAssembler()
    cja.add_statements([st_gef, st_gap])
    cja.make_model()
    assert(len(cja._nodes) == 3)
    assert(len(cja._edges) == 2)
    polarities = [edge['data']['polarity'] for edge in cja._edges]
    assert(len(set(polarities))==2)
    assert('positive' in polarities)
    assert('negative' in polarities)

def test_selfmod():
    cja = CyJSAssembler()
    cja.add_statements([st_selfmod])
    cja.make_model()
    assert(len(cja._nodes) == 1)
    assert(len(cja._edges) == 1)
    polarities = [edge['data']['polarity'] for edge in cja._edges]
    assert(len(polarities) == 1)
    assert(polarities[0] == 'positive')

def test_complex():
    cja = CyJSAssembler()
    cja.add_statements([st_complex])
    cja.make_model()
    assert(len(cja._nodes) == 3)
    assert(len(cja._edges) == 3)
    polarities = [edge['data']['polarity'] for edge in cja._edges]
    assert(len(set(polarities))==1)
    assert('none' in polarities)

def test_print_cyjs_graph():
    cja = CyJSAssembler()
    cja.add_statements([st_act, st_act2])
    cja.make_model()
    cyjs_str = cja.print_cyjs_graph()
    # assert output is not empty
    assert(len(cyjs_str) > len('{\n "edges": [],\n "nodes": []\n}'))

def test_no_grouping():
    st1 = Phosphorylation(Agent('A'), Agent('B'))
    st2 = Phosphorylation(Agent('A'), Agent('C'))
    st3 = Phosphorylation(Agent('C'), Agent('B'))
    cja = CyJSAssembler()
    cja.add_statements([st1, st2, st3])
    cja.make_model(grouping=True)
    parents = [node['data']['parent'] for node in cja._nodes]
    for parent in parents:
        assert parent == ''

def test_grouping_block_targeting_node():
    st1 = Phosphorylation(Agent('A'), Agent('B'))
    st2 = Phosphorylation(Agent('C'), Agent('B'))
    cja = CyJSAssembler()
    cja.add_statements([st1, st2])
    cja.make_model(grouping=True)
    for node in cja._nodes:
        if node['data']['name'] == 'A':
            parent_a = node['data']['parent']
        if node['data']['name'] == 'B':
            parent_b = node['data']['parent']
            assert(parent_b == '')
        if node['data']['name'] == 'C':
            parent_c = node['data']['parent']
    assert_element_properties(cja)
    assert(parent_a == parent_c)
    parent_a_name = [x['data']['name'] for x in cja._nodes if
                     x['data']['id']==parent_a][0]
    assert(parent_a_name.startswith('Group'))
    assert(len(cja._edges) == 3)
    virtual_edges = [x for x in cja._edges if
                     x['data']['i'] == 'Virtual']
    assert(len(virtual_edges) == 2)
    real_edges = [x for x in cja._edges if
                  x['data']['i'] != 'Virtual']
    assert(len(real_edges) == 1)

def test_grouping_node_targeting_block():
    st1 = Phosphorylation(Agent('A'), Agent('B'))
    st2 = Phosphorylation(Agent('A'), Agent('C'))
    cja = CyJSAssembler()
    cja.add_statements([st1, st2])
    cja.make_model(grouping=True)
    for node in cja._nodes:
        if node['data']['name'] == 'A':
            parent_a = node['data']['parent']
            assert(parent_a == '')
        if node['data']['name'] == 'B':
            parent_b = node['data']['parent']
        if node['data']['name'] == 'C':
            parent_c = node['data']['parent']
    assert_element_properties(cja)
    assert(parent_b == parent_c)
    parent_b_name = [x['data']['name'] for x in cja._nodes if
                     x['data']['id']==parent_b][0]
    assert(parent_b_name.startswith('Group'))
    assert(len(cja._edges) == 3)
    virtual_edges = [x for x in cja._edges if
                     x['data']['i'] == 'Virtual']
    assert(len(virtual_edges) == 2)
    real_edges = [x for x in cja._edges if
                  x['data']['i'] != 'Virtual']
    assert(len(real_edges) == 1)

def test_grouping_node_targeting_block_targeting_node():
    st1 = Phosphorylation(Agent('A'), Agent('B'))
    st2 = Phosphorylation(Agent('A'), Agent('C'))
    st3 = Phosphorylation(Agent('B'), Agent('D'))
    st4 = Phosphorylation(Agent('C'), Agent('D'))
    cja = CyJSAssembler()
    cja.add_statements([st1, st2, st3, st4])
    cja.make_model(grouping=True)
    for node in cja._nodes:
        if node['data']['name'] == 'A':
            parent_a = node['data']['parent']
            assert(parent_a == '')
        if node['data']['name'] == 'B':
            parent_b = node['data']['parent']
        if node['data']['name'] == 'C':
            parent_c = node['data']['parent']
        if node['data']['name'] == 'D':
            parent_d = node['data']['parent']
            assert(parent_d == '')
    assert_element_properties(cja)
    assert(parent_b == parent_c)
    parent_b_name = [x['data']['name'] for x in cja._nodes if
                     x['data']['id']==parent_b][0]
    assert(parent_b_name.startswith('Group'))
    assert(len(cja._edges) == 6)
    virtual_edges = [x for x in cja._edges if
                     x['data']['i'] == 'Virtual']
    assert(len(virtual_edges) == 4)
    real_edges = [x for x in cja._edges if
                  x['data']['i'] != 'Virtual']
    assert(len(real_edges) == 2)

def test_grouping_block_targeting_block():
    st1 = Phosphorylation(Agent('A'), Agent('B'))
    st2 = Phosphorylation(Agent('A'), Agent('C'))
    st3 = Phosphorylation(Agent('D'), Agent('B'))
    st4 = Phosphorylation(Agent('D'), Agent('C'))
    cja = CyJSAssembler()
    cja.add_statements([st1, st2, st3, st4])
    cja.make_model(grouping=True)
    for node in cja._nodes:
        if node['data']['name'] == 'A':
            parent_a = node['data']['parent']
        if node['data']['name'] == 'B':
            parent_b = node['data']['parent']
        if node['data']['name'] == 'C':
            parent_c = node['data']['parent']
        if node['data']['name'] == 'D':
            parent_d = node['data']['parent']
    assert_element_properties(cja)
    assert(parent_b == parent_c)
    assert(parent_a == parent_d)
    parent_b_name = [x['data']['name'] for x in cja._nodes if
                     x['data']['id']==parent_b][0]
    parent_a_name = [x['data']['name'] for x in cja._nodes if
                     x['data']['id']==parent_a][0]
    assert(parent_b_name.startswith('Group'))
    assert(parent_a_name.startswith('Group'))
    assert(len(cja._edges) == 5)
    virtual_edges = [x for x in cja._edges if
                     x['data']['i'] == 'Virtual']
    assert(len(virtual_edges) == 4)
    real_edges = [x for x in cja._edges if
                  x['data']['i'] != 'Virtual']
    assert(len(real_edges) == 1)

def test_edge_aggregation_between_nongroup_nodes():
    cja = CyJSAssembler()
    cja.add_statements([st_phos_Y, st_phos_T])
    cja.make_model(grouping=False)
    assert(len(cja._nodes) == 2)
    assert(len(cja._edges) == 1)
    for edge in cja._edges:
        assert(len(edge['data']['uuid_list']) == 2)
    for node in cja._nodes:
        assert(len(node['data']['uuid_list']) == 2)
    cja = CyJSAssembler()
    cja.add_statements([st_phos_Y, st_phos_T])
    cja.make_model(grouping=True)
    assert(len(cja._nodes) == 2)
    assert(len(cja._edges) == 1)
    for edge in cja._edges:
        assert(len(edge['data']['uuid_list']) == 2)
    for node in cja._nodes:
        assert(len(node['data']['uuid_list']) == 2)

def assert_element_properties(cja):
    # each element needs an id
    elements = ([n for n in cja._nodes] + [e for e in cja._edges])
    for element in elements:
        assert(element['data']['id'] is not None), "Element ID is none"
        assert(element['data']['id'] != ''), "Element ID is blank string!"
        # each element should also have a list of uuids with at least one uuid
        assert(element['data']['uuid_list'] is not None), "uuid_list is None"
        assert(len(element['data']['uuid_list']) >= 1), "uuid_list is empty!"
        for uuid in element['data']['uuid_list']:
            assert(type(uuid) == type('abc')), (str(uuid) + ' is not a string')
