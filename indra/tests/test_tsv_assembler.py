import os
from indra.sources import signor
from indra.assemblers.tsv_assembler import TsvAssembler

# Get some statements from Signor
from .test_signor import test_data_file, test_complexes_file
sp = signor.process_from_file(test_data_file, test_complexes_file)
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


