from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import indra.assemblers.english_assembler as ea
from indra.statements import *

def test_agent_basic():
    s = ea._assemble_agent_str(Agent('EGFR'))
    print(s)
    assert (s == 'EGFR')

def test_agent_mod():
    mc = ModCondition('phosphorylation')
    a = Agent('EGFR', mods=mc)
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'phosphorylated EGFR')

def test_agent_mod2():
    mc = ModCondition('phosphorylation', 'tyrosine')
    a = Agent('EGFR', mods=mc)
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'tyrosine-phosphorylated EGFR')

def test_agent_mod3():
    mc = ModCondition('phosphorylation', 'tyrosine', '1111')
    a = Agent('EGFR', mods=mc)
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR phosphorylated on Y1111')

def test_agent_mods():
    mc1 = ModCondition('phosphorylation', 'tyrosine', '1111')
    mc2 = ModCondition('phosphorylation', 'tyrosine', '1234')
    a = Agent('EGFR', mods=[mc1, mc2])
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR phosphorylated on Y1111 and Y1234')

def test_agent_mods2():
    mc1 = ModCondition('phosphorylation', 'tyrosine', '1111')
    mc2 = ModCondition('phosphorylation', 'tyrosine')
    a = Agent('EGFR', mods=[mc1, mc2])
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR phosphorylated on Y1111 and tyrosine')

def test_agent_mods3():
    mc1 = ModCondition('phosphorylation', 'tyrosine', '1111')
    mc2 = ModCondition('phosphorylation')
    a = Agent('EGFR', mods=[mc1, mc2])
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR phosphorylated on Y1111 and an unknown residue')

def test_agent_bound():
    bc = BoundCondition(Agent('EGF'), True)
    a = Agent('EGFR', bound_conditions=[bc])
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR bound to EGF')

def test_agent_not_bound():
    bc = BoundCondition(Agent('EGF'), False)
    a = Agent('EGFR', bound_conditions=[bc])
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR not bound to EGF')

def test_agent_bound_two():
    bc = BoundCondition(Agent('EGF'), True)
    bc2 = BoundCondition(Agent('EGFR'), True)
    a = Agent('EGFR', bound_conditions=[bc, bc2])
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR bound to EGF and EGFR')

def test_agent_bound_three():
    bc = BoundCondition(Agent('EGF'), True)
    bc2 = BoundCondition(Agent('EGFR'), True)
    bc3 = BoundCondition(Agent('GRB2'), True)
    a = Agent('EGFR', bound_conditions=[bc, bc2, bc3])
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR bound to EGF, EGFR and GRB2')

def test_agent_bound_mixed():
    bc = BoundCondition(Agent('EGF'), True)
    bc2 = BoundCondition(Agent('EGFR'), False)
    a = Agent('EGFR', bound_conditions=[bc, bc2])
    s = ea._assemble_agent_str(a)
    print(s)
    assert (s == 'EGFR bound to EGF and not bound to EGFR')

def test_phos_noenz():
    a = Agent('MAP2K1')
    st = Phosphorylation(None, a)
    s = ea._assemble_modification(st)
    print(s)
    assert(s == 'MAP2K1 is phosphorylated.')

def test_phos_noenz2():
    a = Agent('MAP2K1')
    st = Phosphorylation(None, a, 'serine')
    s = ea._assemble_modification(st)
    print(s)
    assert(s == 'MAP2K1 is phosphorylated on serine.')

def test_phos_noenz3():
    a = Agent('MAP2K1')
    st = Phosphorylation(None, a, 'serine', '222')
    s = ea._assemble_modification(st)
    print(s)
    assert(s == 'MAP2K1 is phosphorylated on S222.')

def test_phos_enz():
    a = Agent('MAP2K1')
    b = Agent('BRAF')
    st = Phosphorylation(b, a, 'serine', '222')
    s = ea._assemble_modification(st)
    print(s)
    assert(s == 'BRAF phosphorylates MAP2K1 on S222.')

def test_phos_enz():
    a = Agent('MAP2K1')
    b = Agent('PP2A')
    st = Dephosphorylation(b, a, 'serine', '222')
    s = ea._assemble_modification(st)
    print(s)
    assert(s == 'PP2A dephosphorylates MAP2K1 on S222.')

def test_ubiq_stmt():
    st = Ubiquitination(Agent('X'), Agent('Y'))
    s = ea._assemble_modification(st)
    print(s)
    assert(s == 'X ubiquitinates Y.')

def test_deubiq_stmt():
    st = Deubiquitination(Agent('X'), Agent('Y'))
    s = ea._assemble_modification(st)
    print(s)
    assert(s == 'X deubiquitinates Y.')

def test_deubiq_noenz():
    st = Deubiquitination(None, Agent('Y'))
    s = ea._assemble_modification(st)
    print(s)
    assert(s == 'Y is deubiquitinated.')

def test_complex_one():
    a = Agent('MAP2K1')
    b = Agent('BRAF')
    st = Complex([a, b])
    s = ea._assemble_complex(st)
    print(s)
    assert(s == 'MAP2K1 binds BRAF.')

def test_complex_more():
    a = Agent('MAP2K1')
    b = Agent('BRAF')
    c = Agent('RAF1')
    st = Complex([a, b, c])
    s = ea._assemble_complex(st)
    print(s)
    assert(s == 'MAP2K1 binds BRAF and RAF1.')

def test_assemble_one():
    a = Agent('MAP2K1')
    b = Agent('PP2A')
    st = Dephosphorylation(b, a, 'serine', 222)
    e = ea.EnglishAssembler()
    e.add_statements([st])
    s = e.make_model()
    print(s)
    assert(s == 'PP2A dephosphorylates MAP2K1 on S222.')

def test_assemble_more():
    a = Agent('MAP2K1')
    b = Agent('PP2A')
    st1 = Dephosphorylation(b, a, 'serine', 222)
    b = Agent('BRAF')
    c = Agent('RAF1')
    st2 = Complex([a, b, c])
    e = ea.EnglishAssembler()
    e.add_statements([st1, st2])
    s = e.make_model()
    print(s)
    assert(s ==\
        'PP2A dephosphorylates MAP2K1 on S222. MAP2K1 binds BRAF and RAF1.')

def test_autophos():
    a = Agent('EGFR')
    st = Autophosphorylation(a, 'Y')
    e = ea.EnglishAssembler()
    e.add_statements([st])
    s = e.make_model()
    print(s)
    assert(s == 'EGFR phosphorylates itself on tyrosine.')

def test_activation():
    st = Activation(Agent('MEK'), 'activity',
                    Agent('ERK'), 'activity', True)
    e = ea.EnglishAssembler()
    e.add_statements([st])
    s = e.make_model()
    print(s)
    assert(s == 'MEK activates ERK.')

def test_agent_loc():
    a = Agent('BRAF', location='cytoplasm')
    print(ea._assemble_agent_str(a))
    assert(ea._assemble_agent_str(a) == 'BRAF in the cytoplasm')

def test_translocation():
    st1 = Translocation(Agent('FOXO3A'))
    st2 = Translocation(Agent('FOXO3A'), 'cytoplasm')
    st3 = Translocation(Agent('FOXO3A'), None, 'nucleus')
    st4 = Translocation(Agent('FOXO3A'), 'cytoplasm', 'nucleus')
    e = ea.EnglishAssembler()
    e.add_statements([st1])
    s = e.make_model()
    assert(s == 'FOXO3A translocates.')
    e = ea.EnglishAssembler()
    e.add_statements([st2])
    s = e.make_model()
    assert(s == 'FOXO3A translocates from the cytoplasm.')
    e = ea.EnglishAssembler()
    e.add_statements([st3])
    s = e.make_model()
    assert(s == 'FOXO3A translocates to the nucleus.')
    e = ea.EnglishAssembler()
    e.add_statements([st4])
    s = e.make_model()
    assert(s == 'FOXO3A translocates from the cytoplasm to the nucleus.')
