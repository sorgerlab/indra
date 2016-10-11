from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers import CyJSAssembler

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

def test_complex():
    cja = CyJSAssembler()
    cja.add_statements([st_complex])
    cja.make_model()
    assert(len(cja._nodes) == 3)
    assert(len(cja._edges) == 3)
    polarities = [edge['data']['polarity'] for edge in cja._edges]
    assert(len(set(polarities))==1)
    assert('none' in polarities)

def test_print_cyjs():
    cja = CyJSAssembler()
    cja.add_statements([st_act, st_act2])
    cja.make_model()
    cyjs_str = cja.print_cyjs()
    print(cyjs_str)
