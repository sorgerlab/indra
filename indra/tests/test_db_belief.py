from indra.db.belief import MockStatement, MockEvidence
from indra.belief import BeliefEngine, default_probs


def test_belief_calc():
    be = BeliefEngine()
    test_stmts = [
        MockStatement(
            evidence=[MockEvidence('sparser'), MockEvidence('reach')],
            mk_hash=1,
            supports=[2, 3, 5]
            ),
        MockStatement(
            evidence=MockEvidence('biopax'),
            mk_hash=2,
            supports=[3]
            ),
        MockStatement(
            evidence=MockEvidence('signor'),
            mk_hash=3,
            supports=[]
            ),
        MockStatement(
            evidence=MockEvidence('biogrid'),
            mk_hash=4,
            supports=[3]
            ),
        MockStatement(
            evidence=MockEvidence('bel'),
            mk_hash=5,
            supports=[]
            ),
        MockStatement(
            evidence=[MockEvidence('phosphosite'), MockEvidence('trips')],
            mk_hash=6,
            supports=[]
            )
        ]
    be.set_prior_probs(test_stmts)
    be.set_hierarchy_probs(test_stmts)
    be.set_linked_probs(test_stmts)
    results = {s.matches_key(): s.belief for s in test_stmts}
    print(results)