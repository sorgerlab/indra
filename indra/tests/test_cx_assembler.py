from indra.statements import *
from indra.assemblers.cx import CxAssembler, NiceCxAssembler, hub_layout

mek = Agent('MAP2K1', db_refs={'HGNC': '6840'})
erk = Agent('MAPK1', db_refs={'UP': 'P28482'})
dusp = Agent('DUSP4')
st_phos = Phosphorylation(mek, erk)
st_dephos = Dephosphorylation(dusp, erk)
st_complex = Complex([mek, erk, dusp])
st_complex2 = Complex([mek, mek, erk, erk, dusp])
st_act = Activation(mek, erk)
st_gef = Gef(Agent('SOS1'), Agent('HRAS'))
st_gap = Gap(Agent('RASA1'), Agent('HRAS'))
st_act2 = Inhibition(dusp, erk)
st_not_cited = Phosphorylation(mek, erk, evidence=[Evidence()])
st_cited = Phosphorylation(mek, erk, evidence=[Evidence(pmid='12345',
                                              text='MEK phosphorylates ERK')])
st_invalid_cited = Phosphorylation(mek, erk, evidence=[Evidence(pmid='api35',
                                              text='MEK phosphorylates ERK')])


def test_phos():
    cxa = CxAssembler()
    cxa.add_statements([st_phos])
    cxa.make_model()
    assert len(cxa.cx['nodes']) == 2
    assert len(cxa.cx['edges']) == 1


def test_phos_nice():
    cxa = NiceCxAssembler([st_cited])
    cxa.make_model()
    assert len(cxa.network.nodes) == 2
    assert len(cxa.network.edges) == 1
    print(cxa.print_model())


def test_gapgef_nice():
    cxa = NiceCxAssembler([st_gef, st_gap])
    cxa.make_model()
    assert len(cxa.network.nodes) == 3, cxa.network.nodes
    assert len(cxa.network.edges) == 2
    print(cxa.print_model())


def test_dephos():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    assert len(cxa.cx['nodes']) == 3
    assert len(cxa.cx['edges']) == 2


def test_complex():
    cxa = CxAssembler()
    cxa.add_statements([st_complex])
    cxa.make_model()
    assert len(cxa.cx['nodes']) == 3
    assert len(cxa.cx['edges']) == 3


def test_complex2():
    cxa = CxAssembler()
    cxa.add_statements([st_complex2])
    cxa.make_model()
    assert len(cxa.cx['nodes']) == 3
    assert len(cxa.cx['edges']) == 5


def test_act():
    cxa = CxAssembler()
    cxa.add_statements([st_act, st_act2])
    cxa.make_model()
    assert len(cxa.cx['nodes']) == 3
    assert len(cxa.cx['edges']) == 2


def test_gef():
    cxa = CxAssembler()
    cxa.add_statements([st_gef])
    cxa.make_model()
    assert len(cxa.cx['nodes']) == 2
    assert len(cxa.cx['edges']) == 1


def test_gap():
    cxa = CxAssembler()
    cxa.add_statements([st_gap])
    cxa.make_model()
    assert len(cxa.cx['nodes']) == 2
    assert len(cxa.cx['edges']) == 1


def test_node_attributes():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    assert len(cxa.cx['nodeAttributes']) == 5


def test_edge_attributes():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    assert len(cxa.cx['edgeAttributes']) == 14


def test_cited():
    cxa = CxAssembler()
    cxa.add_statements([st_cited])
    cxa.make_model()
    assert len(cxa.cx['citations']) == 1
    assert len(cxa.cx['edgeCitations']) == 1
    citation = cxa.cx['citations'][0]
    assert citation.get('dc:identifier') == 'pmid:12345'
    cid = citation.get('@id')
    assert cxa.cx['edgeCitations'][0]['citations'][0] == cid
    print(cxa.print_cx())


def test_invalid_cited():
    cxa = CxAssembler()
    cxa.add_statements([st_invalid_cited])
    cxa.make_model()
    assert not cxa.cx['citations']
    assert not cxa.cx['edgeCitations']


def test_supports():
    cxa = CxAssembler()
    cxa.add_statements([st_cited])
    cxa.make_model()
    assert len(cxa.cx['supports']) == 1
    assert len(cxa.cx['edgeSupports']) == 1


def test_set_context():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    cxa.set_context('BT20_BREAST')
    print(cxa.cx['nodeAttributes'])
    assert len(cxa.cx['nodeAttributes']) == 11


def test_make_print_model():
    cxa = CxAssembler()
    cxa.add_statements([st_phos])
    cx_str = cxa.make_model()
    assert cx_str


def test_no_pmid():
    cxa = CxAssembler([st_not_cited])
    cxa.make_model()
    assert not cxa.cx['edgeCitations']


def test_hub_layout():
    stmts = [st_phos, st_dephos, st_act]
    cxa = CxAssembler(stmts)
    cxa.make_model()
    graph = hub_layout.cx_to_networkx(cxa.cx)
    erk = hub_layout.get_node_by_name(graph, 'MAPK1')
    hub_layout.add_semantic_hub_layout(cxa.cx, 'MAPK1')
    assert cxa.cx['cartesianLayout']
    for node in cxa.cx['cartesianLayout']:
        if node['node'] == erk:
            assert node['x'] == 0.0
            assert node['y'] == 0.0
        else:
            assert node['x'] != 0
            assert node['y'] != 0

    node_classes = hub_layout.classify_nodes(graph, erk)
    assert node_classes[hub_layout.get_node_by_name(graph, 'DUSP4')] == \
        (1, 'modification', 'other')
    assert node_classes[hub_layout.get_node_by_name(graph, 'MAP2K1')] in \
        {(1, 'activity', 'protein'), (1, 'modification', 'protein')}
