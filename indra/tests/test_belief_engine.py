from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.belief import BeliefEngine
from indra.belief import _get_belief_package

ev1 = Evidence(source_api='reach')
ev2 = Evidence(source_api='trips')
ev3 = Evidence(source_api='assertion')


def test_prior_prob_one():
    be = BeliefEngine()
    prob = 1 - (be.prior_probs['rand']['reach'] +
                be.prior_probs['syst']['reach'])
    st = Phosphorylation(None, Agent('a'), evidence=[ev1])
    assert(st.belief == 1)
    be.set_prior_probs([st])
    assert(st.belief == prob)


def test_prior_prob_two_same():
    be = BeliefEngine()
    prob = 1 - (be.prior_probs['rand']['reach']**2 +
                be.prior_probs['syst']['reach'])
    st = Phosphorylation(None, Agent('a'), evidence=[ev1, ev1])
    assert(st.belief == 1)
    be.set_prior_probs([st])
    assert(st.belief == prob)


def test_prior_prob_two_different():
    be = BeliefEngine()
    prob = 1 - (be.prior_probs['rand']['reach'] +
                 be.prior_probs['syst']['reach']) * \
               (be.prior_probs['rand']['trips'] +
                 be.prior_probs['syst']['trips'])
    st = Phosphorylation(None, Agent('a'), evidence=[ev1, ev2])
    assert(st.belief == 1)
    be.set_prior_probs([st])
    assert(st.belief == prob)


def test_prior_prob_one_two():
    be = BeliefEngine()
    prob = 1 - (be.prior_probs['rand']['reach']**2 +
                 be.prior_probs['syst']['reach']) * \
               (be.prior_probs['rand']['trips'] +
                 be.prior_probs['syst']['trips'])
    st = Phosphorylation(None, Agent('a'), evidence=[ev1, ev1, ev2])
    assert(st.belief == 1)
    be.set_prior_probs([st])
    assert(st.belief == prob)


def test_prior_prob_assertion():
    be = BeliefEngine()
    st = Phosphorylation(None, Agent('a'), evidence=[ev1, ev1, ev2, ev3])
    assert(st.belief == 1)
    be.set_prior_probs([st])
    assert(st.belief == 1)


def test_hierarchy_probs1():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('b'), evidence=[ev2])
    st2.supports = [st1]
    st1.supported_by = [st2]
    st1.belief = 0.5
    st2.belief = 0.8
    be.set_hierarchy_probs([st1, st2])
    assert(st1.belief == 0.5)
    assert(st2.belief == 0.9)


def test_hierarchy_probs2():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('b'), evidence=[ev2])
    st3 = Phosphorylation(None, Agent('c'), evidence=[ev3])
    st2.supports = [st1]
    st3.supports = [st1, st2]
    st1.supported_by = [st2, st3]
    st2.supported_by = [st3]
    st1.belief = 0.5
    st2.belief = 0.8
    st3.belief = 0.2
    be.set_hierarchy_probs([st1, st2, st3])
    assert(st1.belief == 0.5)
    assert(st2.belief == 0.9)
    assert(st3.belief == 0.92)


def test_hierarchy_probs3():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('b'), evidence=[ev2])
    st3 = Phosphorylation(None, Agent('c'), evidence=[ev3])
    st3.supports = [st1, st2]
    st1.supported_by = [st3]
    st2.supported_by = [st3]
    st1.belief = 0.5
    st2.belief = 0.8
    st3.belief = 0.2
    be.set_hierarchy_probs([st1, st2, st3])
    assert(st1.belief == 0.5)
    assert(st2.belief == 0.8)
    assert(st3.belief == 0.92)


def test_hierarchy_probs4():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('b'), evidence=[ev2])
    st3 = Phosphorylation(None, Agent('c'), evidence=[ev3])
    st4 = Phosphorylation(None, Agent('d'), evidence=[ev1])
    st4.supports = [st1, st2, st3]
    st3.supports = [st1]
    st2.supports = [st1]
    st1.supported_by = [st2, st3, st4]
    st2.supported_by = [st4]
    st3.supported_by = [st4]
    st1.belief = 0.5
    st2.belief = 0.8
    st3.belief = 0.2
    st4.belief = 0.6
    be.set_hierarchy_probs([st1, st2, st3])
    assert(st1.belief == 0.5)
    assert(st2.belief == 0.9)
    assert(st3.belief == 0.6)
    assert(st4.belief == 0.968)


def test_get_belief_package1():
    st1 = Phosphorylation(None, Agent('a'))
    st1.belief = 0.53
    package = _get_belief_package(st1)
    assert(len(package) == 1)
    assert(package[0][0] == 0.53)
    assert(package[0][1] == st1.matches_key())


def test_get_belief_package2():
    st1 = Phosphorylation(None, Agent('A1'))
    st2 = Phosphorylation(None, Agent('A'))
    st1.supported_by = [st2]
    st2.supports = [st1]
    st1.belief = 0.8
    st2.belief = 0.6
    package = _get_belief_package(st1)
    assert(len(package) == 1)
    assert(package[0][0] == 0.8)
    assert(package[0][1] == st1.matches_key())
    package = _get_belief_package(st2)
    assert(len(package) == 2)
    assert(package[0][0] == 0.8)
    assert(package[0][1] == st1.matches_key())
    assert(package[1][0] == 0.6)
    assert(package[1][1] == st2.matches_key())


def test_get_belief_package3():
    st1 = Phosphorylation(Agent('B'), Agent('A1'))
    st2 = Phosphorylation(None, Agent('A1'))
    st3 = Phosphorylation(None, Agent('A'))
    st1.supported_by = [st2, st3]
    st2.supported_by = [st3]
    st2.supports = [st1]
    st3.supports = [st1, st2]
    st1.belief = 0.8
    st2.belief = 0.6
    st3.belief = 0.7
    package = _get_belief_package(st1)
    assert(len(package) == 1)
    assert(set([p[0] for p in package]) == set([0.8]))
    package = _get_belief_package(st2)
    assert(len(package) == 2)
    assert(set([p[0] for p in package]) == set([0.6, 0.8]))
    package = _get_belief_package(st3)
    assert(len(package) == 3)
    assert(set([p[0] for p in package]) == set([0.6, 0.7, 0.8]))

