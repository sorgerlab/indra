from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers.figaro import FigaroAssembler
from indra.statements import Concept, Event, Influence


def test_assemble_influence():
    stmt = Influence(Event(Concept('rainfall')),
                     Event(Concept('crop_yields')))
    fa = FigaroAssembler([stmt])
    fa.make_model()
    assert fa.BN is not None
    assert len(fa.BN.nodes()) == 2
    assert len(fa.BN.edges()) == 1


def test_print_model():
    stmt1 = Influence(Event(Concept('rainfall')),
                      Event(Concept('crop_yields')))
    stmt2 = Influence(Event(Concept('irrigation')),
                      Event(Concept('crop_yields')))
    stmt3 = Influence(Event(Concept('temperature')),
                      Event(Concept('crop_yields')))
    stmt4 = Influence(Event(Concept('rainfall')),
                      Event(Concept('temperature')))
    stmts = [stmt1, stmt2, stmt3, stmt4]
    fa = FigaroAssembler(stmts)
    fa.make_model()
    txt = fa.print_model()
    assert txt is not None
