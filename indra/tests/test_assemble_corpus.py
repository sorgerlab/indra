from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.tools import assemble_corpus as ac
from indra.statements import Activation, Phosphorylation, Agent, Evidence

a = Agent('a', db_refs={'HGNC': '1234', 'TEXT': 'a'})
b = Agent('b', db_refs={'UP': 'P15056', 'TEXT': 'b'})
c = Agent('c', db_refs={'BE': 'XXX', 'TEXT': 'c'})
d = Agent('d', db_refs={'TEXT': 'd'})
e = Agent('e', db_refs={'CHEBI': 'CHEBI:1234', 'TEXT': 'e'})
f = Agent('b', db_refs={'UP': 'P28028', 'TEXT': 'b'})
g = Agent('g', db_refs={'BE': 'ERK'})
h = Agent('g', mods=['x', 'y'], mutations=['x', 'y'], activity='x',
               location='nucleus', bound_conditions=['x', 'y', 'z'])
st1 = Phosphorylation(a, b)
st2 = Phosphorylation(a, d)
st3 = Phosphorylation(c, d)
st4 = Phosphorylation(b, e)
st5 = Phosphorylation(None, b)
st6 = Phosphorylation(None, d)
st7 = Phosphorylation(None, e)
st8 = Phosphorylation(b, f)
st9 = Phosphorylation(None, f)
st10 = Phosphorylation(None, g)
st11 = Phosphorylation(None, h)
st12 = Phosphorylation(a, b, evidence=[Evidence(epistemics={'direct': True})])
st13 = Phosphorylation(a, b, evidence=[Evidence(epistemics={'direct': False})])
st14 = Activation(a, b, 'activity')
st15 = Activation(a, b, 'kinase')
st1.belief = 0.9
st2.belief = 0.8
st3.belief = 0.7

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
    st_out = ac.filter_genes_only([st3], specific_only=True)
    assert len(st_out) == 0

def test_filter_human_only():
    st_out = ac.filter_human_only([st1, st5])
    assert len(st_out) == 2
    st_out = ac.filter_human_only([st8, st9])
    assert len(st_out) == 0

def test_filter_gene_list_one():
    st_out = ac.filter_gene_list([st1, st2], ['a'], 'one')
    assert(len(st_out) == 2)
    st_out = ac.filter_gene_list([st1, st2], ['a'], 'all')
    assert(len(st_out) == 0)
    st_out = ac.filter_gene_list([st1, st2], ['a', 'b'], 'all')
    assert(len(st_out) == 1)

def test_run_preassembly():
    st_out = ac.run_preassembly([st1, st3, st5, st6])
    assert(len(st_out) == 2)

def test_expand_families():
    st_out = ac.expand_families([st10])
    assert(len(st_out) == 2)

def test_strip_agent_context():
    st_out = ac.strip_agent_context([st11])
    assert(len(st_out) == 1)
    assert(not st_out[0].sub.mods)
    assert(not st_out[0].sub.mutations)
    assert(not st_out[0].sub.bound_conditions)
    assert(not st_out[0].sub.activity)
    assert(not st_out[0].sub.location)

def test_filter_direct():
    st_out = ac.filter_direct([st12])
    assert(len(st_out) == 1)
    st_out = ac.filter_direct([st13])
    assert(len(st_out) == 0)

def test_filter_belief():
    st_out = ac.filter_belief([st1, st2, st3], 0.75)
    assert(len(st_out) == 2)

def test_reduce_activities():
    st_out = ac.reduce_activities([st14, st15])
    assert(st_out[0].obj_activity == 'kinase')
    assert(st_out[1].obj_activity == 'kinase')

def test_filter_source():
    ev1 = Evidence(source_api='bel')
    ev2 = Evidence(source_api='biopax')
    ev3 = Evidence(source_api='reach')
    st1 = Activation(Agent('a'), Agent('b'), evidence=[ev3])
    st2  = Activation(Agent('a'), Agent('b'), evidence=[ev1, ev2])
    st3 = Activation(Agent('a'), Agent('b'), evidence=[ev1, ev3])
    st_out = ac.filter_evidence_source([st1, st2], ['reach'], 'one')
    assert(len(st_out) == 1)
    st_out = ac.filter_evidence_source([st1, st2, st3], ['reach'], 'all')
    assert (len(st_out) == 2)
    st_out = ac.filter_evidence_source([st1, st2, st3], ['bel', 'biopax'], 'one')
    assert (len(st_out) == 2)
    st_out = ac.filter_evidence_source([st1, st2, st3], ['bel', 'biopax'], 'all')
    assert (len(st_out) == 1)

def test_map_grounding():
    a = Agent('MEK', db_refs={'TEXT': 'MEK'})
    b = Agent('X', db_refs={'TEXT': 'ERK'})
    st = Activation(a, b)
    st_out = ac.map_grounding([st], do_rename=False)
    assert(len(st_out) == 1)
    assert(st_out[0].subj.db_refs.get('BE'))
    assert(st_out[0].obj.db_refs.get('BE'))
    assert(st_out[0].obj.name == 'X')
    st_out = ac.map_grounding([st], do_rename=True)
    assert(len(st_out) == 1)
    assert(st_out[0].subj.db_refs.get('BE'))
    assert(st_out[0].obj.db_refs.get('BE'))
    assert(st_out[0].obj.name == 'ERK')

def test_map_sequence():
    a = Agent('MAPK1', db_refs={'UP': 'P28482', 'HGNC': '6871'})
    st1 = Phosphorylation(None, a, 'T', '182')
    st2 = Phosphorylation(None, a, 'T', '185')
    st3 = Phosphorylation(None, a, 'Y', '999')
    st_out = ac.map_sequence([st1])
    assert(len(st_out) == 1)
    assert(st_out[0].position == '185')
    st_out = ac.map_sequence([st2])
    assert(len(st_out) == 1)
    assert(st_out[0].position == '185')
    st_out = ac.map_sequence([st3])
    assert(len(st_out) == 0)

def test_filter_by_type():
    st_out = ac.filter_by_type([st1, st14], Phosphorylation)
    assert(len(st_out) == 1)
