from indra.statements import *
from pysb import *

def test_simple_phosphorylation():
    # Create the statement
    a = Agent('A')
    b = Agent('B')
    st = Phosphorylation(a, b, 'T', '185')
    # Now create the PySB model
    Model('test_simple_phosphorylation')
    Monomer('A')
    Monomer('B', ['T185'], {'T185':['u', 'p']})
    Rule('A_phos_B', A() + B(T185='u') >> A() + B(T185='p'),
         Parameter('k', 1))
    Observable('B_phosphoT185', B(T185='p'))
    Initial(A(), Parameter('A_0', 100))
    Initial(B(T185='u'), Parameter('B_0', 100))
    # Now check the model
    mc = ModelChecker(model, [st])
    results = mc.check_model()
    assert len(results) == 1
    assert isinstance(results[0], tuple)
    assert results[0][0] == st
    assert results[0][1] == True

if __name__ == '__main__':
    test_simple_phosphorylation()

