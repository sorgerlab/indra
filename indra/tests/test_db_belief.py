from nose.plugins.attrib import attr

from indra.db.belief import MockStatement, MockEvidence, populate_support, \
    load_mock_statements, calculate_belief
from indra.belief import BeliefEngine
from indra.tests.test_db_client import _PrePaDatabaseTestSetup


def test_belief_calc_up_to_prior():
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


def test_belief_calc_up_to_hierarchy():
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
    init_results = {s.matches_key(): s.belief for s in test_stmts}
    print(init_results)
    supp_links = [(1,2), (1,3), (2,3), (1,5), (4,3)]
    populate_support(test_stmts, supp_links)
    be.set_hierarchy_probs(test_stmts)
    results = {s.matches_key(): s.belief for s in test_stmts}
    print(results)

    # Test a couple very simple properties.
    assert len(results) == len(test_stmts), (len(results), len(test_stmts))
    assert all([0 < b < 1 for b in results.values()]), 'Beliefs out of range.'

    # Test the change from the initial.
    all_deltas_correct = True
    deltas_dict = {}
    for s in test_stmts:
        h = s.matches_key()
        b = s.belief

        # Get results
        res = {'actual': b - init_results[h]}

        # Define expectations.
        if s.supports:
            res['expected'] = 'increase'
            if res['actual'] <= 0:
                all_deltas_correct = False
        else:
            res['expected'] = 'no change'
            if res['actual'] != 0:
                all_deltas_correct = False

        deltas_dict[h] = res
    assert all_deltas_correct, deltas_dict


def _get_prepped_db(num_stmts):
    dts = _PrePaDatabaseTestSetup(num_stmts)
    dts.load_background()
    dts.add_statements()
    dts.insert_pa_statements()
    return dts.test_db


@attr('nonpublic')
def test_mock_stmt_load():
    db = _get_prepped_db(1000)
    stmts = load_mock_statements(db)
    assert 500 <= len(stmts) <= 1000, len(stmts)
    assert all([len(s.evidence) >= 1 for s in stmts])
    sid_list = [ev.annotations['raw_sid'] for s in stmts for ev in s.evidence]
    sid_set = set(sid_list)
    assert len(sid_list) == len(sid_set), (len(sid_list), len(sid_set))
    assert len([sup for s in stmts for sup in s.supports]) \
        == db.count(db.PASupportLinks), "Support is missing."
    belief_dict = calculate_belief(stmts)
    assert len(belief_dict) == len(stmts), (len(belief_dict), len(stmts))
    assert all([0 < b < 1 for b in belief_dict.values()]),\
        'Belief values out of range.'
