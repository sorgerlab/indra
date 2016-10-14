from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.mechlinker import MechLinker

def test_act_phos_to_af():
    act_st = Activation(Agent('A'), 'activity', Agent('B'), 'activity', True)
    phos_st = Phosphorylation(Agent('A'), Agent('B'))
    ml = MechLinker([act_st, phos_st])
    linked_stmts = ml.link_statements()
    assert(len(linked_stmts) == 1)
