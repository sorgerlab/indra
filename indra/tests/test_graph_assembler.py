from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers import GraphAssembler

def test_phosphorylation():
    st = [Phosphorylation(Agent('MAP2K1'), Agent('MAPK1'))]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    assert(len(ga.graph.nodes()) == 2)
    assert(len(ga.graph.edges()) == 1)

def test_phosphorylation_noenz():
    st = [Phosphorylation(None, Agent('MAPK1'))]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    assert(len(ga.graph.nodes()) == 0)
    assert(len(ga.graph.edges()) == 0)

def test_dephosphorylation():
    st = [Dephosphorylation(Agent('DUSP4'), Agent('MAPK1'))]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    assert(len(ga.graph.nodes()) == 2)
    assert(len(ga.graph.edges()) == 1)

def test_dephosphorylation_noenz():
    st = [Dephosphorylation(None, Agent('MAPK1'))]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    assert(len(ga.graph.nodes()) == 0)
    assert(len(ga.graph.edges()) == 0)

def test_activation():
    st = [Activation(Agent('MAP2K1'), 'activity',
                     Agent('MAPK1'), 'activity', True)]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    assert(len(ga.graph.nodes()) == 2)
    assert(len(ga.graph.edges()) == 1)

def test_inactivation():
    st = [Activation(Agent('DUSP4'), 'activity',
                     Agent('MAPK1'), 'activity', False)]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    assert(len(ga.graph.nodes()) == 2)
    assert(len(ga.graph.edges()) == 1)

def test_complex():
    st = [Complex([Agent('BRAF'), Agent('RAF1'), Agent('YWAH')])]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    assert(len(ga.graph.nodes()) == 3)
    assert(len(ga.graph.edges()) == 3)

def test_duplicates():
    st = [Complex([Agent('BRAF'), Agent('RAF1'), Agent('YWAH')])]
    st += [Complex([Agent('BRAF'), Agent('RAF1')])]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    assert(len(ga.graph.nodes()) == 3)
    assert(len(ga.graph.edges()) == 3)

def test_get_string():
    st = [Phosphorylation(Agent('MAP2K1'), Agent('MAPK1'))]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    graph_str = ga.get_string()
    assert(graph_str)

def test_save_dot():
    st = [Phosphorylation(Agent('MAP2K1'), Agent('MAPK1'))]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    ga.save_dot('/dev/null')

def test_save_pdf():
    st = [Phosphorylation(Agent('MAP2K1'), Agent('MAPK1'))]
    ga = GraphAssembler()
    ga.add_statements(st)
    ga.make_model()
    ga.save_pdf('/dev/null')
