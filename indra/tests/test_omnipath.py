from __future__ import unicode_literals
from builtins import dict, str
from indra.statements import Phosphorylation
from indra.databases import omnipath as op

def test_query_ptms():
    stmts = op.get_ptms(['Q13873'])
    assert len(stmts) == 1
    assert isinstance(stmts[0], Phosphorylation)
    assert stmts[0].enz.name == 'CSNK2A1'
    assert stmts[0].sub.name == 'BMPR2'
    assert stmts[0].residue == 'S'
    assert stmts[0].position == '757'
