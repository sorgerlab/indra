from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.util import unicode_strs

ev = Evidence(source_api='bel', pmid='12345', epistemics={'direct': True},
              text='This is the evidence.')

def test_mod_condition_from():
    jd = {'mod_type': 'phosphorylation', 'residue': 'S'}
    mc = ModCondition._from_json(jd)
    assert(mc.residue == 'S')
    assert(mc.mod_type == 'phosphorylation')
    assert(mc.position is None)

def test_agent_mod_condition():
    a = Agent('MAP2K1', mods=[ModCondition('phosphorylation', 'serine', 218),
                              ModCondition('phosphorylation', 'serine', 222)])
    jd = a.to_json()
    jd2 = Agent._from_json(jd).to_json()
    assert(jd == jd2)

def test_modification():
    stmt = Phosphorylation(Agent('a'), Agent('b'), 'S', evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_selfmodification():
    stmt = Autophosphorylation(Agent('a'), 'Y', '1234', evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_activation():
    stmt = Activation(Agent('a'), Agent('b'), 'kinase', evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_amount():
    stmt = IncreaseAmount(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_active_form():
    stmt = ActiveForm(Agent('a', location='nucleus'), 'kinase', False,
                      evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_complex():
    stmt = Complex([Agent('a'), Agent('b')], evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_translocation():
    stmt = Translocation(Agent('a'), 'cytoplasm', 'nucleus', evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_rasgap():
    stmt = RasGap(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_rasgef():
    stmt = RasGef(Agent('a'), Agent('b'), evidence=[ev])
    jd = stmt.to_json()
    jd2 = Statement._from_json(jd).to_json()
    assert(jd == jd2)

def test_supports():
    stmt1 = RasGap(Agent('B'), Agent('B'), evidence=[ev])
    stmt2 = RasGap(Agent('a'), Agent('b'), evidence=[ev])
    stmt1.supports = [stmt2]
    stmt2.supported_by = [stmt1]
    jd1 = stmt1.to_json()
    jd2 = stmt2.to_json()
    jds = [jd1, jd2]
    stmts = stmts_from_json(jds)
    assert(len(stmts[0].supports) == 1)
    assert(len(stmts[1].supported_by) == 1)
    assert(stmts[0].supports[0] == stmts[1])
    assert(stmts[1].supported_by[0] == stmts[0])
    jds2 = stmts_to_json(stmts)
    stmts2 = stmts_from_json(jds2)
    assert(len(stmts2[0].supports) == 1)
    assert(len(stmts2[1].supported_by) == 1)
    assert(stmts2[0].supports[0] == stmts2[1])
    assert(stmts2[1].supported_by[0] == stmts2[0])
