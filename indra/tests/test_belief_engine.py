from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.belief import BeliefEngine

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

'''
def test_hierarchy_probs1():
    be = BeliefEngine()
    st1 = Phosphorylation(None, Agent('a'), evidence=[ev1])
    st2 = Phosphorylation(None, Agent('a'), evidence=[ev2], supports=[st1])
    assert(st.belief == 1)
    be.set_prior_probs([st])
    assert(st.belief == 1)
'''
