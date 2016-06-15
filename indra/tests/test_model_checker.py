from indra.statements import *
from pysb import *
from pysb.core import SelfExporter
from pysb.tools import render_reactions
from indra.tools.model_checker import ModelChecker, mp_embeds_into, \
                                      cp_embeds_into
from pysb.tools import species_graph
from pysb.bng import generate_equations
from pysb import kappa
from pysb.testing import with_model

@with_model
def test_mp_embedding():
    # Create a PySB model
    Monomer('A', ['b', 'other'], {'other':['u','p']})
    mp1 = A(other='u')
    mp2 = A()
    mp3 = A(other='p')
    assert mp_embeds_into(mp1, mp2)
    assert not mp_embeds_into(mp2, mp1)
    assert mp_embeds_into(mp3, mp2)
    assert not mp_embeds_into(mp2, mp3)
    assert not mp_embeds_into(mp3, mp1)
    assert not mp_embeds_into(mp1, mp3)

@with_model
def test_cp_embedding():
    Monomer('A', ['b', 'other'], {'other':['u','p']})
    Monomer('B', ['b'])
    cp1 = A(b=1, other='p') % B(b=1)
    cp2 = A()
    cp3 = A(b=1, other='u') % B(b=1)
    cp4 = A(other='p')
    cp5 = A(b=1) % B(b=1)
    # Some tests not performed because ComplexPatterns for second term are not
    # yet supported
    assert cp_embeds_into(cp1, cp2)
    #assert not cp_embeds_into(cp1, cp3)
    assert cp_embeds_into(cp1, cp4)
    #assert not cp_embeds_into(cp1, cp5)
    #assert not cp_embeds_into(cp2, cp1)
    #assert not cp_embeds_into(cp2, cp3)
    assert not cp_embeds_into(cp2, cp4)
    #assert not cp_embeds_into(cp2, cp5)
    #assert not cp_embeds_into(cp3, cp1)
    assert cp_embeds_into(cp3, cp2)
    assert not cp_embeds_into(cp3, cp4)
    #assert cp_embeds_into(cp3, cp5)
    #assert not cp_embeds_into(cp4, cp1)
    assert cp_embeds_into(cp4, cp2)
    #assert not cp_embeds_into(cp4, cp3)
    #assert not cp_embeds_into(cp4, cp5)
    #assert not cp_embeds_into(cp5, cp1)
    assert cp_embeds_into(cp5, cp2)
    #assert not cp_embeds_into(cp5, cp3)
    assert not cp_embeds_into(cp5, cp4)

"""
@with_model
def test_match_rules():
    Monomer('A')
    Monomer('B', ['T185'], {'T185':['u', 'p']})
    rule = Rule('A_phos_B', A() + B(T185='u') >> A() + B(T185='p'),
                Parameter('k', 1))
"""

@with_model
def test_one_step_phosphorylation():
    # Override the shutoff of self export in psyb_assembler
    # Create the statement
    a = Agent('A')
    b = Agent('B')
    st = Phosphorylation(a, b, 'T', '185')
    # Now create the PySB model
    Monomer('A')
    Monomer('B', ['T185'], {'T185':['u', 'p']})
    Rule('A_phos_B', A() + B(T185='u') >> A() + B(T185='p'),
         Parameter('k', 1))
    Initial(A(), Parameter('A_0', 100))
    Initial(B(T185='u'), Parameter('B_0', 100))
    #with open('model_rxn.dot', 'w') as f:
    #    f.write(render_reactions.run(model))
    #with open('species_1step.dot', 'w') as f:
    #    f.write(species_graph.run(model))
    # Now check the model
    mc = ModelChecker(model, [st])
    results = mc.check_model()
    assert len(results) == 1
    assert isinstance(results[0], tuple)
    assert results[0][0] == st
    assert results[0][1] == True

@with_model
def test_two_step_phosphorylation():
    # Create the statement
    a = Agent('A')
    b = Agent('B')
    st = Phosphorylation(a, b, 'T', '185')
    # Now create the PySB model
    Monomer('A', ['b', 'other'], {'other':['u','p']})
    Monomer('B', ['b', 'T185'], {'T185':['u', 'p']})
    Rule('A_bind_B', A(b=None) + B(b=None, T185='u') <>
                     A(b=1) % B(b=1, T185='u'),
                 Parameter('kf', 1), Parameter('kr', 1))
    Rule('A_phos_B', A(b=1) % B(b=1, T185='u') >>
                     A(b=None) + B(b=None, T185='p'),
                 Parameter('kcat', 1))
    Initial(A(b=None, other='p'), Parameter('Ap_0', 100))
    Initial(A(b=None, other='u'), Parameter('Au_0', 100))
    Initial(B(b=None, T185='u'), Parameter('B_0', 100))
    #with open('model_rxn.dot', 'w') as f:
    #    f.write(render_reactions.run(model))
    #with open('species_2step.dot', 'w') as f:
    #    f.write(species_graph.run(model))
    #im = kappa.influence_map(model)
    #im.draw('im_2step.pdf', prog='dot')
    #generate_equations(model)
    # Now check the model
    mc = ModelChecker(model, [st])
    results = mc.check_model()
    assert len(results) == 1
    assert isinstance(results[0], tuple)
    assert results[0][0] == st
    assert results[0][1] == True


if __name__ == '__main__':
    test_mp_embedding()
    test_cp_embedding()

