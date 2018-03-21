import os
from indra.assemblers.tsv_assembler import TsvAssembler
from indra.sources.signor import SignorProcessor

# Get some statements from Signor
#sp = SignorProcessor()
#stmts = sp.statements

from indra.tools import assemble_corpus as ac
#ac.dump_statements(stmts, 'test_stmts.pkl')
stmts = ac.load_statements('test_stmts.pkl')

def test_tsv_init():
    ta = TsvAssembler(stmts)
    ta.make_model('tsv_test')

def test_tsv_add_stmts():
    ta = TsvAssembler()
    ta.add_statements(stmts)
    assert len(ta.statements) == len(stmts)

def test_make_model():
    ta = TsvAssembler(stmts)
    ta.make_model('tsv_test')
    assert os.path.exists('tsv_test_mods.tsv')

if __name__ == '__main__':
    test_tsv_init()
    test_tsv_add_stmts()
    test_make_model()


