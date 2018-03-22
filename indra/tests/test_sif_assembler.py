from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers import SifAssembler
from indra.statements import *


def test_simple_assembly():
    st1 = Activation(Agent('a'), Agent('b'))
    st2 = Inhibition(Agent('a'), Agent('c'))
    sa = SifAssembler([st1, st2])
    sa.make_model()
    assert(len(sa.graph.nodes()) == 3)
    assert(len(sa.graph.edges()) == 2)


def test_evidence_assembly():
    ev1 = Evidence(pmid='1')
    ev2 = Evidence(pmid='2')
    ev3 = Evidence(pmid='3')
    ev4 = Evidence(pmid='4')
    st1 = Activation(Agent('a'), Agent('b'), evidence=[ev1])
    st2 = Inhibition(Agent('a'), Agent('c'), evidence=[ev1, ev2, ev3])
    sa = SifAssembler([st1, st2])
    sa.make_model()
    assert(len(sa.graph.nodes()) == 3)
    assert(len(sa.graph.edges()) == 2)
    sa.set_edge_weights('support_pmid')


def test_modification():
    st1 = Phosphorylation(Agent('BRAF'), Agent('MAP2K1'), 'S', '222')
    sa = SifAssembler([st1])
    sa.make_model(True, True, True)
    assert(len(sa.graph.nodes()) == 2)
    assert(len(sa.graph.edges()) == 1)
    sa.save_model('test_sif.sif', True)
    with open('test_sif.sif', 'rb') as fh:
        txt = fh.read().decode('utf-8')
        assert txt == 'BRAF 0 MAP2K1\n', txt

