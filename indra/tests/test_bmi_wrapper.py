from indra.statements import *
from indra.assemblers import PysbAssembler
from indra.assemblers.bmi_wrapper import BMIModel

stmts = [Influence(Concept('rainfall'), Concept('flood')),
         Influence(Concept('flood'), Concept('displacement'))]


def make_bmi_model():
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    bm = BMIModel(model)
    return bm


def test_bmi_model():
    bm = make_bmi_model()
    assert len(bm.model.monomers) == 3
    assert len(bm.model.rules) == 2


def test_initialize():
    bm = make_bmi_model()
    bm.initialize()
    assert bm.sim is not None
    assert bm.time == 0.0


def test_update():
    bm = make_bmi_model()
    bm.initialize()
    bm.update(dt=100)
    assert bm.time == 100.0
