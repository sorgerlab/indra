import pickle
from indra.tools import assemble_corpus as ac
from indra.statements import Phosphorylation, Agent

a = Agent('a', db_refs={'HGNC': '1234', 'TEXT': 'a'})
b = Agent('b', db_refs={'UP': 'P15056', 'TEXT': 'b'})
c = Agent('c', db_refs={'BE': 'XXX', 'TEXT': 'c'})
d = Agent('d', db_refs={'TEXT': 'd'})
e = Agent('e', db_refs={'CHEBI': 'CHEBI:1234', 'TEXT': 'e'})
f = Agent('b', db_refs={'UP': 'P28028', 'TEXT': 'b'})
st1 = Phosphorylation(a, b)
st2 = Phosphorylation(a, d)
st3 = Phosphorylation(c, d)
st4 = Phosphorylation(b, e)
st5 = Phosphorylation(None, b)
st6 = Phosphorylation(None, d)
st7 = Phosphorylation(None, e)
st8 = Phosphorylation(b, f)
st9 = Phosphorylation(None, f)

def test_load_stmts():
    with open('_test.pkl', 'wb') as fh:
        pickle.dump([st1], fh, protocol=2)
    st_loaded = ac.load_statements('_test.pkl')
    assert(len(st_loaded) == 1)
    assert(st_loaded[0].equals(st1))

def test_dump_stmts():
    ac.dump_statements([st1], '_test.pkl')
    st_loaded = ac.load_statements('_test.pkl')
    assert(len(st_loaded) == 1)
    assert(st_loaded[0].equals(st1))

def test_filter_grounded_only():
    st_out = ac.filter_grounded_only([st1, st4])
    assert len(st_out) == 2
    st_out = ac.filter_grounded_only([st3])
    assert len(st_out) == 0

def test_filter_genes_only():
    st_out = ac.filter_genes_only([st1, st5])
    assert len(st_out) == 2
    st_out = ac.filter_genes_only([st6, st7])
    assert len(st_out) == 0
    st_out = ac.filter_genes_only([st4])
    assert len(st_out) == 0

def test_filter_human_only():
    st_out = ac.filter_human_only([st1, st5])
    assert len(st_out) == 2
    st_out = ac.filter_human_only([st8, st9])
    assert len(st_out) == 0

def test_run_preassembly():
    st_out = ac.run_preassembly([st1, st3, st5, st6])
    assert(len(st_out) == 2)
