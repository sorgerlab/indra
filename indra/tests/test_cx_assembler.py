from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers import CxAssembler

mek = Agent('MAP2K1', db_refs={'HGNC': '6840'})
erk = Agent('MAPK1', db_refs={'UP': 'P28482'})
dusp = Agent('DUSP4')
st_phos = Phosphorylation(mek, erk)
st_dephos = Dephosphorylation(dusp, erk)
st_complex = Complex([mek, erk, dusp])
st_act = Activation(mek, 'activity', erk, 'activity', True)
st_rasgef = RasGef(Agent('SOS1'), 'activity', Agent('HRAS'))
st_rasgap = RasGap(Agent('RASA1'), 'activity', Agent('HRAS'))
st_act2 = Activation(dusp, 'activity', erk, 'activity', False)
st_cited = Phosphorylation(mek, erk, evidence=Evidence(pmid='12345',
                                              text='MEK phosphorylates ERK'))
st_cited2 = Phosphorylation(mek, erk, evidence=Evidence(pmid='api35',
                                              text='MEK phosphorylates ERK'))

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

def test_act():
    cxa = CxAssembler()
    cxa.add_statements([st_act, st_act2])
    cxa.make_model()
    assert(len(cxa.cx['nodes']) == 3)
    assert(len(cxa.cx['edges']) == 2)

def test_rasgef():
    cxa = CxAssembler()
    cxa.add_statements([st_rasgef])
    cxa.make_model()
    assert(len(cxa.cx['nodes']) == 2)
    assert(len(cxa.cx['edges']) == 1)

def test_rasgap():
    cxa = CxAssembler()
    cxa.add_statements([st_rasgap])
    cxa.make_model()
    assert(len(cxa.cx['nodes']) == 2)
    assert(len(cxa.cx['edges']) == 1)

def test_node_attributes():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    assert(len(cxa.cx['nodeAttributes']) == 5)

def test_edge_attributes():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    assert(len(cxa.cx['edgeAttributes']) == 8)

def test_cited():
    cxa = CxAssembler()
    cxa.add_statements([st_cited])
    cxa.make_model()
    assert(len(cxa.cx['citations']) == 1)
    assert(len(cxa.cx['edgeCitations']) == 1)

def test_cited():
    cxa = CxAssembler()
    cxa.add_statements([st_cited2])
    cxa.make_model()
    assert(len(cxa.cx['citations']) == 1)
    assert(len(cxa.cx['edgeCitations']) == 1)

def test_supports():
    cxa = CxAssembler()
    cxa.add_statements([st_cited])
    cxa.make_model()
    assert(len(cxa.cx['supports']) == 1)
    assert(len(cxa.cx['edgeSupports']) == 1)

def test_set_context():
    cxa = CxAssembler()
    cxa.add_statements([st_phos, st_dephos])
    cxa.make_model()
    cxa.set_context('BT20_BREAST')
    assert(len(cxa.cx['nodeAttributes']) == 8)
