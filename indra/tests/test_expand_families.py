from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
from indra.preassembler.hierarchy_manager import hierarchies

from indra.statements import Agent, Phosphorylation
from indra.tools import expand_families as ef


def test_expand_families():
    raf = Agent('RAF', db_refs={'BE':'RAF'})
    mek = Agent('MEK', db_refs={'BE':'MEK'})
    st = Phosphorylation(raf, mek)
    exp = ef.Expander(hierarchies)
    expanded_stmts = exp.expand_families([st])

