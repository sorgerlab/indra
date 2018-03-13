from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers import CAGAssembler

eg1 = {'EIDOS': [('a/b/c', 0.123)]}
eg2 = {'EIDOS': [('a/b/c', 0.234)]}

statements = [Influence(
    Agent('inorganic fertilizer', db_refs=eg1),
    Agent('farm sizes', db_refs=eg2),
    {'adjectives': 'serious', 'polarity': 1},
    {'adjectives': 'significant', 'polarity': 1},
)]


def test_assemble_influence():
    ca = CAGAssembler(statements)
    CAG = ca.make_model()
    assert(len(CAG.nodes()) == 2)
    assert(len(CAG.edges()) == 1)


def test_export_to_cyjs():
    ca = CAGAssembler(statements)
    ca.make_model()
    cyjs = ca.export_to_cytoscapejs()
    assert len(cyjs['nodes']) == 2
    assert len(cyjs['edges']) == 1
