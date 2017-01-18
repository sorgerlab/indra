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
