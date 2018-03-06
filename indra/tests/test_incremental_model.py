from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.tools.incremental_model import IncrementalModel

stmts = [Complex([Agent('MAPK1'), Agent('MAPK3')]),
         Complex([Agent('MAPK3'), Agent('MAP2K1')])]
stmts2 = [Phosphorylation(None, Agent('MAPK3', db_refs={'UP': '123'}))]
stmt3 = Phosphorylation(None, Agent('BRAF', db_refs={'HGNC': '1097',
                                                     'UP': 'P15056'}))
stmt4 = Phosphorylation(None, Agent('RAF', db_refs={'FPLX': 'RAF'}))
stmt5 = Phosphorylation(Agent('X'), Agent('RAF', db_refs={'FPLX': 'RAF'}))

def test_add_stmts_blank():
    im = IncrementalModel()
    im.add_statements('12345', stmts)
    assert(len(im.get_statements()) == 2)
    im.preassemble()
    assert(len(im.assembled_stmts) == 2)

def test_add_stmts_blank_nofilter():
    im = IncrementalModel()
    im.add_statements('12345', stmts)
    im.preassemble(filters=None)
    assert(len(im.assembled_stmts) == 2)

def test_add_stmts_blank_emptyfilter():
    im = IncrementalModel()
    im.add_statements('12345', stmts)
    im.preassemble(filters=[])
    assert(len(im.assembled_stmts) == 2)

def test_add_stmts_blank_noprior():
    im = IncrementalModel()
    im.add_statements('12345', stmts)
    im.preassemble(filters=['prior_one'])
    assert(len(im.assembled_stmts) == 2)

def test_add_stmts_blank_noprior2():
    im = IncrementalModel()
    im.add_statements('12345', stmts)
    im.preassemble(filters=['prior_all'])
    assert(len(im.assembled_stmts) == 2)

def test_add_stmts_prior_one():
    im = IncrementalModel()
    im.stmts['prior'] = [stmts[0]]
    im.prior_genes = ['MAPK1', 'MAPK3']
    im.add_statements('12345', [stmts[1]])
    im.preassemble(filters=['prior_one'])
    assert(len(im.assembled_stmts) == 2)

def test_add_stmts_prior_all():
    im = IncrementalModel()
    im.stmts['prior'] = [stmts[0]]
    im.prior_genes = ['MAPK1', 'MAPK3']
    im.add_statements('12345', [stmts[1]])
    im.preassemble(filters=['prior_all'])
    assert(len(im.assembled_stmts) == 1)

def test_grounding_not_all():
    im = IncrementalModel()
    stmt = Complex([Agent('A', db_refs={'UP': 'ABCD'}), 
                   Agent('B')])
    im.add_statements('12345', [stmt])
    im.preassemble(filters=['grounding'])
    assert(len(im.assembled_stmts) == 0)

def test_grounding_all():
    im = IncrementalModel()
    stmt = Complex([Agent('A', db_refs={'UP': 'ABCD'}), 
                   Agent('B', db_refs={'HGNC': '1234'})])
    im.add_statements('12345', [stmt])
    im.preassemble(filters=['grounding'])
    assert(len(im.assembled_stmts) == 1)

def test_grounding_none():
    im = IncrementalModel()
    im.add_statements('12345', stmts)
    im.preassemble(filters=['grounding'])
    assert(len(im.assembled_stmts) == 0)

def test_grounding_none_agent():
    im = IncrementalModel()
    im.add_statements('12345', stmts2)
    im.preassemble(filters=['grounding'])
    assert(len(im.assembled_stmts) == 1)

def test_human_only():
    im = IncrementalModel()
    stmt1 = Phosphorylation(None, Agent('BRAF', db_refs={'UP': 'P15056'}))
    stmt2 = Phosphorylation(None, Agent('BRAF', db_refs={'UP': 'P28028'}))
    stmt3 = Phosphorylation(None, Agent('BRAF', db_refs={'HGNC': 'BRAF'}))
    stmt4 = Phosphorylation(None, Agent('BRAF', db_refs={}))
    im.add_statements('12345', [stmt1])
    im.preassemble(filters=['human_only'])
    assert(len(im.assembled_stmts) == 1)
    im.add_statements('12346', [stmt2])
    im.preassemble(filters=['human_only'])
    assert(len(im.assembled_stmts) == 1)
    im.add_statements('12346', [stmt3])
    im.preassemble(filters=['human_only'])
    assert (len(im.assembled_stmts) == 2)
    im.add_statements('12346', [stmt4])
    im.preassemble(filters=['human_only'])
    assert (len(im.assembled_stmts) == 3)
