import os
from indra.assemblers.tsv_assembler import TsvAssembler
from indra.sources.signor import SignorProcessor

# Get some statements from Signor
sp = SignorProcessor()
stmts = sp.statements

def test_tsv_init():
    ta = TsvAssembler(stmts)
    ta.make_model('tsv_test')

def test_tsv_add_stmts():
    ta = TsvAssembler()
    ta.add_statements(stmts)
    assert len(ta.statements) == len(stmts)

def test_make_model():
    ta = TsvAssembler(stmts)
    ta.make_model('tsv_test.tsv')
    assert os.path.exists('tsv_test.tsv')


