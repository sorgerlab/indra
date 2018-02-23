from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers import FigaroAssembler
from indra.statements import Agent, Influence


def test_assemble_influence():
    stmt = Influence(Agent('rainfall'), Agent('crop_yields'))
    fa = FigaroAssembler([stmt])
    fa.make_model()
    assert fa.BN is not None
    assert len(fa.BN.nodes()) == 2
    assert len(fa.BN.edges()) == 1


def test_print_model():
    stmt1 = Influence(Agent('rainfall'), Agent('crop_yields'))
    stmt2 = Influence(Agent('irrigation'), Agent('crop_yields'))
    stmt3 = Influence(Agent('temperature'), Agent('crop_yields'))
    stmt4 = Influence(Agent('rainfall'), Agent('temperature'))
    stmts = [stmt1, stmt2, stmt3, stmt4]
    fa = FigaroAssembler(stmts)
    fa.make_model()
    txt = fa.print_model()
    assert txt is not None
