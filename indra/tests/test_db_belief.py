from indra.db.belief import MockStatement, MockEvidence
from indra.belief import BeliefEngine, default_probs


def test_belief_calc_no_supps():
    be = BeliefEngine()
    test_stmts = [
        MockStatement(1, [MockEvidence('sparser'), MockEvidence('reach')]),
        MockStatement(2, MockEvidence('biopax')),
        MockStatement(3, MockEvidence('signor')),
        MockStatement(4, MockEvidence('biogrid')),
        MockStatement(5, MockEvidence('bel')),
        MockStatement(6, [MockEvidence('phosphosite'), MockEvidence('trips')]),
        ]
    be.set_prior_probs(test_stmts)
    results = {s.matches_key(): s.belief for s in test_stmts}
    print(results)
    assert len(results) == len(test_stmts), (len(results), len(test_stmts))
    assert all([0 < b < 1 for b in results.values()]), 'Beliefs out of range.'
