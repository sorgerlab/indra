import json
from indra.statements import *
from indra.cx_assembler import CxAssembler

mek = Agent('MAP2K1', db_refs={'HGNC': '6840'})
erk = Agent('MAPK1', db_refs={'UP': 'P28482'})
dusp = Agent('DUSP4')
st_phos = Phosphorylation(mek, erk)
st_dephos = Dephosphorylation(dusp, erk)
st_complex = Complex([mek, erk, dusp])
st_actact = ActivityActivity(mek, 'Activity', 'increases', erk, 'Activity')
st_actact2 = ActivityActivity(dusp, 'Activity', 'decreases', erk, 'Activity')
st_cited = Phosphorylation(mek, erk, evidence=Evidence(pmid='12345'))


def test_phos():
    cxa = CxAssembler()
    cxa.add_statements([st_phos])
    cxa.make_model()
    assert(len(cxa.cx['nodes']) == 2)
    assert(len(cxa.cx['edges']) == 1)

def test_dephos():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    assert(len(cxa.cx['nodes']) == 3)
    assert(len(cxa.cx['edges']) == 2)

def test_complex():
    cxa = CxAssembler()
    cxa.add_statements([st_complex])
    cxa.make_model()
    assert(len(cxa.cx['nodes']) == 3)
    assert(len(cxa.cx['edges']) == 3)

def test_actact():
    cxa = CxAssembler()
    cxa.add_statements([st_actact, st_actact2])
    cxa.make_model()
    assert(len(cxa.cx['nodes']) == 3)
    assert(len(cxa.cx['edges']) == 2)

def test_node_attributes():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    assert(len(cxa.cx['nodeAttributes']) == 2)

def test_edge_attributes():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    assert(len(cxa.cx['edgeAttributes']) == 2)
    cxa.save_model('small_example.cx')
