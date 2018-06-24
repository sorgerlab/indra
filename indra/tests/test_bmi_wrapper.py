from indra.statements import *
from indra.assemblers import PysbAssembler
from indra.assemblers.bmi_wrapper import BMIModel

stmts = [Influence(Concept('rainfall'), Concept('flood')),
         Influence(Concept('flood'), Concept('displacement'))]


def make_bmi_model():
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    bm = BMIModel(model, inputs=['rainfall'])
    return bm


def test_bmi_model():
    bm = make_bmi_model()
    assert len(bm.model.monomers) == 3
    assert len(bm.model.rules) == 2


def test_initialize():
    bm = make_bmi_model()
    bm.initialize()
    assert bm.state[1] == 10000.0
    assert bm.sim is not None
    assert bm.time == 0.0


def test_get_in_out_vars():
    bm = make_bmi_model()
    bm.initialize()
    assert set(bm.get_output_var_names()) == {'flood', 'displacement'}
    assert bm.get_input_var_names() == ['rainfall']


def test_update():
    bm = make_bmi_model()
    bm.initialize()
    bm.update(dt=100)
    assert bm.time == 100.0
    assert bm.state[1] != 0.0
    assert bm.get_current_time() == 100.0


def test_set_value():
    bm = make_bmi_model()
    bm.initialize()
    bm.set_value('rainfall', 10.0)
    assert bm.state[bm.species_name_map['rainfall']] == 10.0


def test_get_value():
    bm = make_bmi_model()
    bm.initialize()
    bm.set_value('rainfall', 10.0)
    val = bm.get_value('rainfall')
    assert val == 10.0


def test_get_attribute():
    bm = make_bmi_model()
    attr = bm.get_attribute('model_name')
    assert attr == 'indra_model', attr


def test_make_repo_component():
    bm = make_bmi_model()
    bm.model.name = 'indra_model'
    comp = bm.make_repository_component()
    print(comp)
    assert '<class_name>' in comp
