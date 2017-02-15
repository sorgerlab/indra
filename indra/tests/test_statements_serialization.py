from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.util import unicode_strs

ev = Evidence(source_api='bel', pmid='12345', epistemics={'direct': True},
              text='This is the evidence.')

def test_mod_condition_from():
    jd = {'mod_type': 'phosphorylation', 'residue': 'S'}
    mc = ModCondition.from_json(jd)
    assert(mc.residue == 'S')
    assert(mc.mod_type == 'phosphorylation')
    assert(mc.position is None)

def test_agent_mod_condition():
    a = Agent('MAP2K1', mods=[ModCondition('phosphorylation', 'serine', 218),
                              ModCondition('phosphorylation', 'serine', 222)])
    jd = a.to_json()
    jd2 = Agent.from_json(jd).to_json() 
    assert(jd == jd2)

def test_modification():
    stmt = Phosphorylation(Agent('a'), Agent('b'), 'S', evidence=[ev])
    jd = stmt.to_json()
    jd2 = Phosphorylation.from_json(jd).to_json()
    assert(jd == jd2)

def test_activation():
    stmt = Activation(Agent('a'), Agent('b'), 'kinase', evidence=[ev])
    jd = stmt.to_json()
    jd2 = Activation.from_json(jd).to_json()
    assert(jd == jd2)

def test_amount():
    stmt = IncreaseAmount(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    jd2 = IncreaseAmount.from_json(jd).to_json()
    assert(jd == jd2)

def test_active_form():
    stmt = ActiveForm(Agent('a', location='nucleus'), 'kinase', False,
                      evidence=[ev])
    jd = stmt.to_json()
    jd2 = ActiveForm.from_json(jd).to_json()
    assert(jd == jd2)

def test_complex():
    stmt = Complex([Agent('a'), Agent('b')], evidence=[ev])
    jd = stmt.to_json()
    jd2 = Complex.from_json(jd).to_json()
    assert(jd == jd2)
