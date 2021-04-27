from copy import deepcopy
from nose.tools import raises
from indra.statements import *
from indra.belief import BeliefEngine, load_default_probs, \
    sample_statements, evidence_random_noise_prior, tag_evidence_subtype, \
    SimpleScorer
from indra.belief import BayesianScorer

default_probs = load_default_probs()

ev1 = Evidence(source_api='reach')
ev2 = Evidence(source_api='trips')
ev3 = Evidence(source_api='assertion')
ev4 = Evidence(source_api='biopax')


def test_prior_prob_one():
    be = BeliefEngine()
    prob = 1 - (default_probs['rand']['reach'] +
                default_probs['syst']['reach'])
    st = Phosphorylation(None, Agent('a'), evidence=[ev1])
    assert st.belief == 1
    be.set_prior_probs([st])
    assert st.belief == prob


def test_prior_prob_two_same():
    be = BeliefEngine()
    prob = 1 - (default_probs['rand']['reach']**2 +
                default_probs['syst']['reach'])
    st = Phosphorylation(None, Agent('a'), evidence=[ev1, deepcopy(ev1)])
    assert st.belief == 1
    be.set_prior_probs([st])
    assert st.belief == prob


def test_prior_prob_two_different():
    be = BeliefEngine()
    prob = 1 - (default_probs['rand']['reach'] +
                default_probs['syst']['reach']) * \
               (default_probs['rand']['trips'] +
                default_probs['syst']['trips'])
    st = Phosphorylation(None, Agent('a'), evidence=[ev1, ev2])
    assert st.belief == 1
    be.set_prior_probs([st])
    assert st.belief == prob


def test_prior_prob_one_two():
    be = BeliefEngine()
    prob = 1 - (default_probs['rand']['reach']**2 +
                default_probs['syst']['reach']) * \
               (default_probs['rand']['trips'] +
                default_probs['syst']['trips'])
    st = Phosphorylation(None, Agent('a'), evidence=[ev1, deepcopy(ev1), ev2])
    assert st.belief == 1
    be.set_prior_probs([st])
    assert st.belief == prob


def test_prior_prob_assertion():
    be = BeliefEngine()
    st = Phosphorylation(None, Agent('a'),
                         evidence=[ev1, deepcopy(ev1), ev2, ev3])
    assert st.belief == 1
    be.set_prior_probs([st])
    assert st.belief == 1


def test_hierarchy_probs1():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('b'), evidence=[ev2])
    st2.supports = [st1]
    st1.supported_by = [st2]
    be.set_hierarchy_probs([st1, st2])
    assert_close_enough(st1.belief, 1-0.35)
    assert_close_enough(st2.belief, 1-0.35*0.35)


def test_hierarchy_probs2():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('b'), evidence=[ev2])
    st3 = Phosphorylation(None, Agent('c'), evidence=[ev4])
    st2.supports = [st1]
    st3.supports = [st1, st2]
    st1.supported_by = [st2, st3]
    st2.supported_by = [st3]
    be.set_hierarchy_probs([st1, st2, st3])
    assert_close_enough(st1.belief, 1-0.35)
    assert_close_enough(st2.belief, 1-0.35*0.35)
    assert_close_enough(st3.belief, 1-0.35*0.35*0.21)


def test_hierarchy_probs3():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('b'), evidence=[ev2])
    st3 = Phosphorylation(None, Agent('c'), evidence=[ev4])
    st3.supports = [st1, st2]
    st1.supported_by = [st3]
    st2.supported_by = [st3]
    be.set_hierarchy_probs([st1, st2, st3])
    assert_close_enough(st1.belief, 1-0.35)
    assert_close_enough(st2.belief, 1-0.35)
    assert_close_enough(st3.belief, 1-0.35*0.35*0.21)


def test_hierarchy_probs4():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('b'), evidence=[ev2])
    st3 = Phosphorylation(None, Agent('c'), evidence=[deepcopy(ev1)])
    st4 = Phosphorylation(None, Agent('d'), evidence=[deepcopy(ev1)])
    st4.supports = [st1, st2, st3]
    st3.supports = [st1]
    st2.supports = [st1]
    st1.supported_by = [st2, st3, st4]
    st2.supported_by = [st4]
    st3.supported_by = [st4]
    be.set_hierarchy_probs([st1, st2, st3, st4])
    assert_close_enough(st1.belief, 1-0.35)
    assert_close_enough(st2.belief, 1-0.35*0.35)
    assert_close_enough(st3.belief, 1-(0.05 + 0.3*0.3))
    assert_close_enough(st4.belief, 1-0.35*(0.05 + 0.3*0.3*0.3))


def test_default_probs():
    """Make sure default probs are set with empty constructor."""
    be = BeliefEngine()
    for err_type in ('rand', 'syst'):
        for k, v in default_probs[err_type].items():
            assert default_probs[err_type][k] == v


def test_default_probs_override():
    """Make sure default probs are overriden by constructor argument."""
    prior_probs = {'rand': {'assertion': 0.5}}
    scorer = SimpleScorer(prior_probs)

    be = BeliefEngine(scorer)
    for err_type in ('rand', 'syst'):
        for k, v in scorer.prior_probs[err_type].items():
            if err_type == 'rand' and k == 'assertion':
                assert v == 0.5
            else:
                assert default_probs[err_type][k] == v


def test_default_probs_extend():
    """Make sure default probs are extended by constructor argument."""
    prior_probs = {'rand': {'new_source': 0.1},
                   'syst': {'new_source': 0.05}}
    scorer = SimpleScorer(prior_probs)

    be = BeliefEngine(scorer)
    for err_type in ('rand', 'syst'):
        assert 'new_source' in scorer.prior_probs[err_type]
        for k, v in scorer.prior_probs[err_type].items():
            if err_type == 'rand' and k == 'new_source':
                assert v == 0.1
            elif err_type == 'syst' and k == 'new_source':
                assert v == 0.05
            else:
                assert default_probs[err_type][k] == v


def test_sample_statements():
    st1 = Phosphorylation(Agent('B'), Agent('A1'))
    st2 = Phosphorylation(None, Agent('A1'))
    st3 = Phosphorylation(None, Agent('A'))
    st1.supported_by = [st2, st3]
    st2.supported_by = [st3]
    st2.supports = [st1]
    st3.supports = [st1, st2]
    st1.belief = 0.8
    st2.belief = 0.6
    st3.belief = 0.5
    stmts = sample_statements([st1, st2, st3], seed=10)
    # Random numbers are 0.77132064  0.02075195  0.63364823 with this seed
    assert len(stmts) == 2
    assert st1 in stmts
    assert st2 in stmts
    assert st3 not in stmts


@raises(Exception)
def test_check_prior_probs():
    be = BeliefEngine()
    st = Phosphorylation(None, Agent('ERK'),
                         evidence=[Evidence(source_api='xxx')])
    be.set_prior_probs([st])


def test_evidence_subtype_tagger():
    #Test for reach evidence
    evidence_reach = Evidence(source_api='reach', source_id=0,
            pmid=0, text=None, epistemics={},
            annotations={'found_by': 'Positive_early_activation'})
    (stype, subtype) = tag_evidence_subtype(evidence_reach)
    assert stype == 'reach'
    assert subtype == 'Positive_early_[^_]*'

    #Test for biopax evidence
    evidence_biopax = Evidence(source_api='biopax', source_id=0,
            pmid=0, text=None, epistemics={},
            annotations={'source_sub_id': 'reactome'})
    (stype, subtype) = tag_evidence_subtype(evidence_biopax)
    assert stype == 'biopax'
    assert subtype == 'reactome'

    #Test for geneways evidence
    evidence_geneways = Evidence(source_api='geneways', source_id=0,
            pmid=0, text=None, epistemics={},
            annotations={'actiontype': 'bind'})
    (stype, subtype) = tag_evidence_subtype(evidence_geneways)
    assert stype == 'geneways'
    assert subtype == 'bind'

    #Test for unsupported evidence
    evidence_donald_duck = Evidence(source_api='donald_duck', source_id=0,
            pmid=29053813, text=None, epistemics={'direct': True},
            annotations={'quack': 'quack',
                         'quack?' : 'QUAAAAAAAAACK!'})
    (stype, subtype) = tag_evidence_subtype(evidence_donald_duck)
    assert stype == 'donald_duck'
    assert subtype is None


def test_evidence_random_noise_prior():
    type_probs = {'biopax': 0.9, 'geneways': 0.2}
    biopax_subtype_probs = {
            'reactome': 0.4,
            'biogrid': 0.2}
    geneways_subtype_probs = {
            'phosphorylate': 0.5,
            'bind': 0.7}
    subtype_probs = {'biopax': biopax_subtype_probs,
                     'geneways': geneways_subtype_probs}

    ev_geneways_bind = Evidence(source_api='geneways', source_id=0,
                                pmid=0, text=None, epistemics={},
                                annotations={'actiontype': 'bind'})
    ev_biopax_reactome = Evidence(source_api='biopax', source_id=0,
                                  pmid=0, text=None, epistemics={},
                                  annotations={'source_sub_id': 'reactome'})
    ev_biopax_pid = Evidence(source_api='biopax', source_id=0,
                             pmid=0, text=None, epistemics={},
                             annotations={'source_sub_id': 'pid'})

    # Random noise prior for geneways bind evidence is the subtype prior,
    # since we specified it
    assert evidence_random_noise_prior(ev_geneways_bind, \
                                       type_probs, subtype_probs) == 0.7

    # Random noise prior for reactome biopax evidence is the subtype prior,
    # since we specified it
    assert evidence_random_noise_prior(ev_biopax_reactome, \
                                       type_probs, subtype_probs) == 0.4

    # Random noise prior for pid evidence is the subtype prior,
    # since we specified it
    assert evidence_random_noise_prior(ev_biopax_pid,
                                       type_probs, subtype_probs) == 0.9

    # Make sure this all still works when we go through the belief engine
    statements = []
    members = [Agent('a'), Agent('b')]
    statements.append(Complex(members, evidence=ev_geneways_bind))
    statements.append(Complex(members, evidence=ev_biopax_reactome))
    statements.append(Complex(members, evidence=ev_biopax_pid))
    p = {'rand': type_probs, 'syst': {'biopax': 0, 'geneways': 0}}

    scorer = SimpleScorer(p, subtype_probs)
    engine = BeliefEngine(scorer)
    engine.set_prior_probs(statements)
    assert statements[0].belief == 1 - 0.7
    assert statements[1].belief == 1 - 0.4
    assert statements[2].belief == 1 - 0.9


def test_negative_evidence():
    prior_probs = {'rand': {'new_source': 0.1},
                   'syst': {'new_source': 0.05}}
    getev = lambda x: Evidence(source_api='new_source',
                               epistemics={'negated': x})
    evs1 = [getev(x) for x in [True, True, False]]
    evs2 = [getev(x) for x in [False, False, False]]
    evs3 = [getev(x) for x in [True, True, True]]
    stmts = [Phosphorylation(None, Agent('a'), evidence=e)
             for e in [evs1, evs2, evs3]]
    scorer = SimpleScorer(prior_probs)
    engine = BeliefEngine(scorer)
    engine.set_prior_probs(stmts)
    pr = prior_probs['rand']['new_source']
    ps = prior_probs['syst']['new_source']
    assert_close_enough(stmts[0].belief, ((1-pr)-ps)*(1-((1-pr*pr)-ps)))
    assert_close_enough(stmts[1].belief, (1-pr*pr*pr)-ps)
    assert stmts[2].belief == 0


def test_bayesian_scorer():
    prior_counts = {'hume': [3, 1]}
    subtype_counts = {'eidos': {'rule1': [2, 2], 'rule2': [1, 4]}}
    scorer = BayesianScorer(prior_counts, subtype_counts)
    # Check initial probability assignment
    assert scorer.prior_probs['rand']['hume'] == 0.2
    assert scorer.prior_probs['syst']['hume'] == 0.05
    assert scorer.subtype_probs['eidos']['rule1'] == 0.45
    assert scorer.subtype_probs['eidos']['rule2'] == 0.75
    # Now try to do some updates
    scorer.update_counts({'hume': [0, 2]}, {})
    assert scorer.prior_counts['hume'] == [3, 3]
    scorer.update_counts({}, {'eidos': {'rule1': [6, 0]}})
    assert scorer.subtype_counts['eidos']['rule1'] == [8, 2]
    # Now check that the probabilities are up to date
    assert scorer.prior_probs['rand']['hume'] == 0.45
    assert scorer.prior_probs['syst']['hume'] == 0.05
    assert_close_enough(scorer.subtype_probs['eidos']['rule1'],
                        0.15)
    assert scorer.subtype_probs['eidos']['rule2'] == 0.75


@raises(AssertionError)
def test_cycle():
    st1 = Phosphorylation(Agent('B'), Agent('A1'))
    st2 = Phosphorylation(None, Agent('A1'))
    st1.supports = [st2]
    st1.supported_by = [st2]
    st2.supports = [st1]
    st2.supported_by = [st1]
    engine = BeliefEngine()
    engine.set_hierarchy_probs([st1, st2])


def assert_close_enough(b1, b2):
    assert abs(b1 - b2) < 1e-6, 'Got %.6f, Expected: %.6f' % (b1, b2)
