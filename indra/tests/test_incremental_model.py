from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.tools.incremental_model import IncrementalModel

stmts = [Complex([Agent('A'), Agent('B')]), Complex([Agent('B'), Agent('C')])]
stmts2 = [Phosphorylation(None, Agent('B', db_refs={'UP': '123'}))]

def test_add_stmts_blank():
    im = IncrementalModel()
    im.add_statements('12345', stmts)
    assert(len(im.get_statements()) == 2)

def test_add_stmts_blank_nofilter():
    im = IncrementalModel()
    im.add_statements('12345', stmts, filters=None)
    assert(len(im.get_statements()) == 2)

def test_add_stmts_blank_emptyfilter():
    im = IncrementalModel()
    im.add_statements('12345', stmts, filters=[])
    assert(len(im.get_statements()) == 2)

def test_add_stmts_blank_noprior():
    im = IncrementalModel()
    im.add_statements('12345', stmts, filters=['prior_one'])
    assert(len(im.get_statements()) == 2)

def test_add_stmts_blank_noprior2():
    im = IncrementalModel()
    im.add_statements('12345', stmts, filters=['prior_all'])
    assert(len(im.get_statements()) == 2)

def test_add_stmts_model_one():
    im = IncrementalModel()
    im.add_statements('12345', [stmts[0]])
    im.add_statements('23456', [stmts[1]], filters=['model_one'])
    assert(len(im.get_statements()) == 2)

def test_add_stmts_model_all():
    im = IncrementalModel()
    im.add_statements('12345', [stmts[0]])
    im.add_statements('23456', [stmts[1]], filters=['model_all'])
    assert(len(im.get_statements()) == 1)

def test_add_stmts_prior_one():
    im = IncrementalModel()
    im.stmts['prior'] = [stmts[0]]
    im.add_statements('12345', [stmts[1]], filters=['prior_one'])
    assert(len(im.get_statements()) == 2)

def test_add_stmts_prior_all():
    im = IncrementalModel()
    im.stmts['prior'] = [stmts[0]]
    im.add_statements('12345', [stmts[1]], filters=['prior_all'])
    assert(len(im.get_statements()) == 1)

def test_grounding_not_all():
    im = IncrementalModel()
    stmt = Complex([Agent('A', db_refs={'UP': 'ABCD'}), 
                   Agent('B')])
    im.add_statements('12345', [stmt], filters=['grounding'])
    assert(len(im.get_statements()) == 0)

def test_grounding_all():
    im = IncrementalModel()
    stmt = Complex([Agent('A', db_refs={'UP': 'ABCD'}), 
                   Agent('B', db_refs={'HGNC': '1234'})])
    im.add_statements('12345', [stmt], filters=['grounding'])
    assert(len(im.get_statements()) == 1)

def test_grounding_none():
    im = IncrementalModel()
    im.add_statements('12345', stmts, filters=['grounding'])
    assert(len(im.get_statements()) == 0)

def test_grounding_none_agent():
    im = IncrementalModel()
    im.add_statements('12345', stmts2, filters=['grounding'])
    assert(len(im.get_statements()) == 1)

def test_human_only():
    im = IncrementalModel()
    stmt1 = Phosphorylation(None, Agent('BRAF', db_refs={'UP': 'P15056'}))
    stmt2 = Phosphorylation(None, Agent('BRAF', db_refs={'UP': 'P28028'}))
    im.add_statements('12345', [stmt1], filters=['human_only'])
    assert(len(im.get_statements()) == 1)
    im.add_statements('12346', [stmt2], filters=['human_only'])
    assert(len(im.get_statements()) == 1)
