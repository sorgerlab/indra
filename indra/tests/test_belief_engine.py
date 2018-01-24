from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from nose.tools import raises
from indra.statements import *
from indra.belief import BeliefEngine
from indra.belief import _get_belief_package, default_probs, \
        sample_statements, evidence_random_noise_prior, tag_evidence_subtype

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

def test_default_probs():
    """Make sure default probs are set with empty constructor."""
    be = BeliefEngine()
    for err_type in ('rand', 'syst'):
        for k, v in be.prior_probs[err_type].items():
            assert default_probs[err_type][k] == v

def test_default_probs_override():
    """Make sure default probs are overriden by constructor argument."""
    be = BeliefEngine(prior_probs={'rand': {'assertion': 0.5}})
    for err_type in ('rand', 'syst'):
        for k, v in be.prior_probs[err_type].items():
            if err_type == 'rand' and k == 'assertion':
                assert v == 0.5
            else:
                assert default_probs[err_type][k] == v

def test_default_probs_extend():
    """Make sure default probs are extended by constructor argument."""
    be = BeliefEngine(prior_probs={'rand': {'new_source': 0.1},
                                   'syst': {'new_source': 0.05}})
    for err_type in ('rand', 'syst'):
        assert 'new_source' in be.prior_probs[err_type]
        for k, v in be.prior_probs[err_type].items():
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
    assert(stype == 'reach')
    assert(subtype == 'Positive_early_[^_]*')

    #Test for biopax evidence
    evidence_biopax = Evidence(source_api='biopax', source_id=0,
            pmid=0, text=None, epistemics={},
            annotations={'source_sub_id': 'reactome'})
    (stype, subtype) = tag_evidence_subtype(evidence_biopax)
    assert(stype == 'biopax')
    assert(subtype == 'reactome')

    #Test for geneways evidence
    evidence_geneways = Evidence(source_api='geneways', source_id=0,
            pmid=0, text=None, epistemics={},
            annotations={'actiontype': 'bind'})
    (stype, subtype) = tag_evidence_subtype(evidence_geneways)
    assert(stype == 'geneways')
    assert(subtype == 'bind')

    #Test for unsupported evidence
    evidence_donald_duck = Evidence(source_api='donald_duck', source_id=0,
            pmid=29053813, text=None, epistemics={'direct': True},
            annotations={'quack': 'quack',
                         'quack?' : 'QUAAAAAAAAACK!'})
    (stype, subtype) = tag_evidence_subtype(evidence_donald_duck)
    assert(stype == 'donald_duck')
    assert(subtype is None)

def test_evidence_random_noise_prior():
    type_probs = {'biopax': 0.9, 'geneways': 0.2}
    biopax_subtype_probs = {
            'reactome' : 0.4,
            'biogrid' : 0.2}
    geneways_subtype_probs = {
            'phosphorylate': 0.5,
            'bind' : 0.7}
    subtype_probs = {
        'biopax' : biopax_subtype_probs,
        'geneways' : geneways_subtype_probs }

    ev_geneways_bind = Evidence(source_api='geneways', source_id=0,
            pmid=0, text=None, epistemics={},
            annotations={'actiontype': 'bind'})
    ev_biopax_reactome = Evidence(source_api='biopax', source_id=0,
            pmid=0, text=None, epistemics={},
            annotations={'source_sub_id': 'reactome'})
    ev_biopax_pid = Evidence(source_api='biopax', source_id=0,
            pmid=0, text=None, epistemics={},
            annotations={'source_sub_id': 'pid'})

    #Random noise prior for geneways bind evidence is the subtype prior,
    #since we specified it
    assert(evidence_random_noise_prior(ev_geneways_bind,
        type_probs, subtype_probs) == 0.7)

    #Random noise prior for reactome biopax evidence is the subtype prior,
    #since we specified it
    assert(evidence_random_noise_prior(ev_biopax_reactome,
        type_probs, subtype_probs) == 0.4)

    #Random noise prior for pid evidence is the subtype prior,
    #since we specified it
    assert(evidence_random_noise_prior(ev_biopax_pid,
        type_probs, subtype_probs) == 0.9)

    #Make sure this all still works when we go through the belief engine
    statements = []
    members = [Agent('a'), Agent('b')]
    statements.append( Complex(members, evidence=ev_geneways_bind) )
    statements.append( Complex(members, evidence=ev_biopax_reactome) )
    statements.append( Complex(members, evidence=ev_biopax_pid) )
    p = {'rand': type_probs, 'syst': {'biopax':0, 'geneways':0}}
    engine = BeliefEngine(p, subtype_probs)
    engine.set_prior_probs(statements)
    assert(statements[0].belief == 1 - 0.7)
    assert(statements[1].belief == 1 - 0.4)
    assert(statements[2].belief == 1 - 0.9)

