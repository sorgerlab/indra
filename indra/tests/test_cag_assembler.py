from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers import CAGAssembler

statements = [Influence(
    Agent('inorganic fertilizer'),
    Agent('farm sizes'),
    {'adjectives': 'serious', 'polarity': 1},
    {'adjectives': 'significant', 'polarity': 1},
)]


def test_influence():
    ca = CAGAssembler(statements)
    CAG = ca.make_model()
    assert(len(CAG.nodes()) == 2)
    assert(len(CAG.edges()) == 1)
