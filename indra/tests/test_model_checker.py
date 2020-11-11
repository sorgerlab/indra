from __future__ import absolute_import, print_function, unicode_literals

import os
import unittest
from collections import Counter

import numpy as np
import networkx as nx
from indra.statements import *
from pysb import *
from indra.databases import hgnc_client
from indra.explanation.model_checker import PysbModelChecker, \
    UnsignedGraphModelChecker, SignedGraphModelChecker, PybelModelChecker, \
    PathResult
from indra.explanation.model_checker.model_checker import \
    signed_edges_to_signed_nodes
from indra.explanation.model_checker.pysb import _mp_embeds_into, \
    _cp_embeds_into, _match_lhs, remove_im_params
from indra.explanation.reporting import stmt_from_rule, stmts_from_pysb_path, \
    stmts_from_pybel_path, stmts_from_indranet_path, PybelEdge, \
    pybel_edge_to_english
from indra.assemblers.pysb.assembler import PysbAssembler, \
    set_base_initial_condition
from indra.assemblers.indranet import IndraNetAssembler
from indra.assemblers.pybel.assembler import PybelAssembler, _get_agent_node
from pysb.testing import with_model


# Test PysbModelChecker
@with_model
def test_mp_embedding():
    # Create a PySB model
    Monomer('A', ['b', 'other'], {'other':['u','p']})
    mp1 = A(other='u')
    mp2 = A()
    mp3 = A(other='p')
    assert _mp_embeds_into(mp1, mp2)
    assert not _mp_embeds_into(mp2, mp1)
    assert _mp_embeds_into(mp3, mp2)
    assert not _mp_embeds_into(mp2, mp3)
    assert not _mp_embeds_into(mp3, mp1)
    assert not _mp_embeds_into(mp1, mp3)


@with_model
def test_cp_embedding():
    Monomer('A', ['b', 'other'], {'other': ['u','p']})
    Monomer('B', ['b'])
    cp1 = A(b=1, other='p') % B(b=1)
    cp2 = A()
    cp3 = A(b=1, other='u') % B(b=1)
    cp4 = A(other='p')
    cp5 = A(b=1) % B(b=1)
    # FIXME Some tests not performed because ComplexPatterns for second term
    # FIXME are not yet supported
    assert _cp_embeds_into(cp1, cp2)
    #assert not _cp_embeds_into(cp1, cp3)
    assert _cp_embeds_into(cp1, cp4)
    #assert not _cp_embeds_into(cp1, cp5)
    #assert not _cp_embeds_into(cp2, cp1)
    #assert not _cp_embeds_into(cp2, cp3)
    assert not _cp_embeds_into(cp2, cp4)
    #assert not _cp_embeds_into(cp2, cp5)
    #assert not _cp_embeds_into(cp3, cp1)
    assert _cp_embeds_into(cp3, cp2)
    assert not _cp_embeds_into(cp3, cp4)
    #assert _cp_embeds_into(cp3, cp5)
    #assert not _cp_embeds_into(cp4, cp1)
    assert _cp_embeds_into(cp4, cp2)
    #assert not _cp_embeds_into(cp4, cp3)
    #assert not _cp_embeds_into(cp4, cp5)
    #assert not _cp_embeds_into(cp5, cp1)
    assert _cp_embeds_into(cp5, cp2)
    #assert not _cp_embeds_into(cp5, cp3)
    assert not _cp_embeds_into(cp5, cp4)


@with_model
def test__match_lhs():
    Monomer('A', ['other'], {'other': ['u', 'p']})
    Monomer('B', ['T185'], {'T185': ['u', 'p']})
    rule = Rule('A_phos_B', A() + B(T185='u') >> A() + B(T185='p'),
                Parameter('k', 1))
    matching_rules = _match_lhs(A(), model.rules)
    assert len(matching_rules) == 1
    assert matching_rules[0] == rule
    matching_rules = _match_lhs(A(other='u'), model.rules)
    assert len(matching_rules) == 0


"""
@with_model
def test_match_rhs():
    Monomer('A', ['other'], {'other':['u', 'p']})
    Monomer('B', ['T185'], {'T185':['u', 'p']})
    rule = Rule('A_phos_B', A() + B(T185='u') >> A() + B(T185='p'),
                Parameter('k', 1))
    matching_rules = _match_rhs(B(T185='p'), model.rules)
    assert len(matching_rules) == 1
    assert matching_rules[0] == rule
    matching_rules = _match_rhs(B(T185='u'), model.rules)
    assert len(matching_rules) == 0
    matching_rules = _match_rhs(B(), model.rules)
    assert len(matching_rules) == 1
    assert matching_rules[0] == rule
"""


@with_model
def test_one_step_phosphorylation():
    # Create the statement
    a = Agent('A', db_refs={'HGNC': '1'})
    b = Agent('B', db_refs={'HGNC': '2'})
    st = Phosphorylation(a, b, 'T', '185')
    # Now create the PySB model
    Monomer('A')
    Monomer('B', ['T185'], {'T185': ['u', 'p']})
    Rule('A_phos_B', A() + B(T185='u') >> A() + B(T185='p'),
         Parameter('k', 1))
    Initial(A(), Parameter('A_0', 100))
    Initial(B(T185='u'), Parameter('B_0', 100))
    # Add annotations
    Annotation(A, 'https://identifiers.org/hgnc:1')
    Annotation(B, 'https://identifiers.org/hgnc:2')
    Annotation('A_phos_B', 'A', 'rule_has_subject')
    Annotation('A_phos_B', 'B', 'rule_has_object')
    B.site_annotations = [
        Annotation(('T185', 'p'), 'phosphorylation', 'is_modification'),
        Annotation('T185', 'T', 'is_residue'),
        Annotation('T185', '185', 'is_position'),
    ]
    mc = PysbModelChecker(model, [st])
    results = mc.check_model()
    assert len(results) == 1
    assert isinstance(results[0], tuple)
    assert results[0][0] == st
    pr = results[0][1]
    assert isinstance(pr, PathResult)
    assert pr.paths == [(('A_phos_B', 0), ('B_T185_p_obs', 0))]


@with_model
def test_two_step_phosphorylation():
    # Create the statement
    a = Agent('A', db_refs={'HGNC': '1'})
    b = Agent('B', db_refs={'HGNC': '2'})
    st = Phosphorylation(a, b, 'T', '185')
    # Now create the PySB model
    Monomer('A', ['b', 'other'], {'other': ['u','p']})
    Monomer('B', ['b', 'T185'], {'T185': ['u', 'p']})
    Rule('A_bind_B', A(b=None) + B(b=None, T185='u') >>
                     A(b=1) % B(b=1, T185='u'), Parameter('kf', 1))
    Rule('A_bind_B_rev', A(b=1) % B(b=1, T185='u') >>
                         A(b=None) + B(b=None, T185='u'), Parameter('kr', 1))
    Rule('A_phos_B', A(b=1) % B(b=1, T185='u') >>
                     A(b=None) + B(b=None, T185='p'),
                 Parameter('kcat', 1))
    Initial(A(b=None, other='p'), Parameter('Ap_0', 100))
    Initial(A(b=None, other='u'), Parameter('Au_0', 100))
    Initial(B(b=None, T185='u'), Parameter('B_0', 100))
    # Add annotations
    Annotation(A, 'https://identifiers.org/hgnc:1')
    Annotation(B, 'https://identifiers.org/hgnc:2')
    Annotation('A_phos_B', 'A', 'rule_has_subject')
    Annotation('A_phos_B', 'B', 'rule_has_object')
    B.site_annotations = [
        Annotation(('T185', 'p'), 'phosphorylation', 'is_modification'),
        Annotation('T185', 'T', 'is_residue'),
        Annotation('T185', '185', 'is_position'),
    ]
    #with open('model_rxn.dot', 'w') as f:
    #    f.write(render_reactions.run(model))
    #with open('species_2step.dot', 'w') as f:
    #    f.write(species_graph.run(model))
    #generate_equations(model)
    # Now check the model
    mc = PysbModelChecker(model, [st])
    results = mc.check_model()
    assert len(results) == 1
    assert isinstance(results[0], tuple)
    assert results[0][0] == st
    pr = results[0][1]
    assert pr.paths == [(('A_phos_B', 0), ('B_T185_p_obs', 0))]


def test_pysb_assembler_phospho_policies():
    a = Agent('A', db_refs={'HGNC': '1'})
    b = Agent('B', db_refs={'HGNC': '2'})
    st = Phosphorylation(a, b, 'T', '185')
    pa = PysbAssembler()
    pa.add_statements([st])
    # Try two step
    pa.make_model(policies='two_step')
    mc = PysbModelChecker(pa.model, [st])
    results = mc.check_model()
    assert len(results) == 1
    assert isinstance(results[0], tuple)
    assert results[0][0] == st
    pr = results[0][1]
    assert pr.paths == [(('A_phosphorylation_B_T185', 0), ('B_T185_p_obs', 0))]
    assert stmts_from_pysb_path(pr.paths[0], pa.model, [st]) == [st]
    # Try one step
    pa.make_model(policies='one_step')
    mc = PysbModelChecker(pa.model, [st])
    results = mc.check_model()
    assert len(results) == 1
    assert isinstance(results[0], tuple)
    assert results[0][0] == st
    pr = results[0][1]
    assert pr.path_found
    assert pr.paths == [(('A_phosphorylation_B_T185', 0), ('B_T185_p_obs', 0))]
    assert stmts_from_pysb_path(pr.paths[0], pa.model, [st]) == [st]
    # Try interactions_only
    pa.make_model(policies='interactions_only')
    mc = PysbModelChecker(pa.model, [st])
    results = mc.check_model()
    assert len(results) == 1
    assert isinstance(results[0], tuple)
    assert results[0][0] == st
    pr = results[0][1]
    assert not pr.path_found


"""
def test_ras_220_network():
    ras_220_results_path = os.path.join('../../models/ras_220_genes'
                                        '/ras_220_gn_related2_stmts.pkl')
    #ras_220_results_path = 'braf_dusp6_stmts.pkl'
    #ras_220_results_path = 'braf_dusp6_small.pkl'
    with open(ras_220_results_path, 'rb') as f:
        ras220_stmts = pickle.load(f)
    ras220_stmts = [s for s in ras220_stmts
                      if isinstance(s, Modification) or
                         isinstance(s, ActiveForm)]
    print("Done loading")
    # Build a PySB model from the Ras 220 statements
    pa = PysbAssembler()
    pa.add_statements(ras220_stmts)
    pa.make_model(policies='one_step')
    # Now create an indirect statement to check the model against
    egfr = Agent('EGFR')
    braf = Agent('BRAF')
    dusp6 = Agent('DUSP6')
    stmt1 = Phosphorylation(braf, dusp6, 'S', '159')
    stmt2 = Phosphorylation(egfr, dusp6, 'S', '159')
    # Check model
    stmts = [stmt1, stmt2]
    mc = ModelChecker(pa.model, stmts)
    checks = mc.check_model()
    assert len(checks) == 2
    assert isinstance(checks[0], tuple)
    assert checks[0][0] == stmt1
    assert checks[0][1] == True
    assert checks[1][0] == stmt2
    assert checks[1][1] == False
    # Now try again, with a two_step policy
    # Skip this, building the influence map takes a very long time
    #pa.make_model(policies='two_step')
    #mc = ModelChecker(pa.model, [stmt1, stmt2])
    #checks = mc.check_model()
    #print checks
    #assert len(checks) == 2
    #assert isinstance(checks[0], tuple)
    #assert checks[0][0] == stmt1
    #assert checks[0][1] == True
    #assert checks[1][0] == stmt2
    #assert checks[1][1] == False
    # Now with an interactions_only policy
    pa.make_model(policies='interactions_only')
    mc = ModelChecker(pa.model, [stmt1, stmt2])
    checks = mc.check_model()
    assert len(checks) == 2
    assert isinstance(checks[0], tuple)
    assert checks[0][0] == stmt1
    assert checks[0][1] == False
    assert checks[1][0] == stmt2
    assert checks[1][1] == False
"""

"""
def test_path_polarity():
    im = pgv.AGraph('im_polarity.dot')
    path1 = ['BRAF_phospho_MAPK1_T185_1', 'MAPK1_phospho_DUSP6_S159_1']
    path2 = ['BRAF_phospho_MAPK1_T185_1', 'BRAF_phospho_MAPK1_T185_3',
             'MAPK1_phospho_DUSP6_S159_1']
    assert _positive_path(im, path1)
    assert not _positive_path(im, path2)
"""


@with_model
def test_consumption_rule():
    pvd = Agent('Pervanadate', db_refs={'HGNC': '1'})
    erk = Agent('MAPK1', db_refs={'HGNC': '2'})
    stmt = Phosphorylation(pvd, erk, 'T', '185')
    # Now make the model
    Monomer('Pervanadate', ['b'])
    Monomer('DUSP', ['b'])
    Monomer('MAPK1', ['b', 'T185'], {'T185': ['u', 'p']})
    Rule('Pvd_binds_DUSP',
         Pervanadate(b=None) + DUSP(b=None) >>
         Pervanadate(b=1) % DUSP(b=1),
         Parameter('k1', 1))
    Rule('Pvd_binds_DUSP_rev',
         Pervanadate(b=1) % DUSP(b=1) >>
         Pervanadate(b=None) + DUSP(b=None),
         Parameter('k2', 1))
    Rule('DUSP_binds_MAPK1_phosT185',
         DUSP(b=None) + MAPK1(b=None, T185='p') >>
         DUSP(b=1) % MAPK1(b=1, T185='p'),
         Parameter('k3', 1))
    Rule('DUSP_binds_MAPK1_phosT185_rev',
         DUSP(b=1) % MAPK1(b=1, T185='p') >>
         DUSP(b=None) + MAPK1(b=None, T185='p'),
         Parameter('k4', 1))
    Rule('DUSP_dephos_MAPK1_at_T185',
         DUSP(b=1) % MAPK1(b=1, T185='p') >>
         DUSP(b=None) % MAPK1(b=None, T185='u'),
         Parameter('k5', 1))
    Annotation(Pervanadate, 'https://identifiers.org/hgnc:1')
    Annotation(MAPK1, 'https://identifiers.org/hgnc:2')
    Annotation('Pvd_binds_DUSP', 'Pervanadate', 'rule_has_subject')
    Annotation('Pvd_binds_DUSP', 'Pervanadate', 'rule_has_object')
    Annotation('Pvd_binds_DUSP', 'DUSP', 'rule_has_subject')
    Annotation('Pvd_binds_DUSP', 'DUSP', 'rule_has_object')
    Annotation('Pvd_binds_DUSP_rev', 'Pervanadate', 'rule_has_subject')
    Annotation('Pvd_binds_DUSP_rev', 'Pervanadate', 'rule_has_object')
    Annotation('Pvd_binds_DUSP_rev', 'DUSP', 'rule_has_subject')
    Annotation('Pvd_binds_DUSP_rev', 'DUSP', 'rule_has_object')
    Annotation('DUSP_dephos_MAPK1_at_T185', 'DUSP', 'rule_has_subject')
    Annotation('DUSP_dephos_MAPK1_at_T185', 'MAPK1', 'rule_has_object')
    MAPK1.site_annotations = [
            Annotation(('T185', 'p'), 'phosphorylation', 'is_modification'),
            Annotation('T185', 'T', 'is_residue'),
            Annotation('T185', '185', 'is_position'),
        ]
    # Now check the model against the statement
    mc = PysbModelChecker(model, [stmt])
    checks = mc.check_model()
    assert len(checks) == 1
    assert isinstance(checks[0], tuple)
    assert checks[0][0] == stmt
    pr = checks[0][1]
    assert pr.paths == [(('Pvd_binds_DUSP', 0),
                         ('DUSP_binds_MAPK1_phosT185', 1),
                         ('DUSP_dephos_MAPK1_at_T185', 1),
                         ('MAPK1_T185_p_obs', 0))]


def test_dephosphorylation():
    dusp = Agent('DUSP6', db_refs={'HGNC':'1'})
    mapk1 = Agent('MAPK1', db_refs={'HGNC':'2'})
    stmt = Dephosphorylation(dusp, mapk1, 'T', '185')

    def check_policy(policy, result):
        pysba = PysbAssembler()
        pysba.add_statements([stmt])
        pysba.make_model(policies=policy)
        mc = PysbModelChecker(pysba.model, [stmt])
        checks = mc.check_model()
        assert len(checks) == 1
        assert isinstance(checks[0], tuple)
        assert checks[0][0] == stmt
        pr = checks[0][1]
        assert pr.paths == result

    check_policy('one_step', [(('DUSP6_dephosphorylation_MAPK1_T185', 0),
                               ('MAPK1_T185_p_obs', 1))])
    check_policy('two_step', [(('DUSP6_dephosphorylation_MAPK1_T185', 0),
                               ('MAPK1_T185_p_obs', 1))])
    check_policy('interactions_only', [])


@with_model
def test_invalid_modification():
    # Override the shutoff of self export in psyb_assembler
    # Create the statement
    a = Agent('A')
    b = Agent('B')
    st = Phosphorylation(a, b, 'T', '185')
    # Now create the PySB model
    Monomer('A')
    Monomer('B', ['Y187'], {'Y187':['u', 'p']})
    Rule('A_phos_B', A() + B(Y187='u') >> A() + B(Y187='p'),
         Parameter('k', 1))
    #Initial(A(), Parameter('A_0', 100))
    #Initial(B(T187='u'), Parameter('B_0', 100))
    #with open('model_rxn.dot', 'w') as f:
    #    f.write(render_reactions.run(model))
    #with open('species_1step.dot', 'w') as f:
    #    f.write(species_graph.run(model))
    # Now check the model
    mc = PysbModelChecker(model, [st])
    results = mc.check_model()
    assert len(results) == 1
    #assert isinstance(results[0], tuple)
    #assert results[0][0] == st
    #assert results[0][1] == True


def _path_polarity_stmt_list():
    a = Agent('A', db_refs={'HGNC': '1'})
    b = Agent('B', db_refs={'HGNC': '2'})
    c = Agent('C', db_refs={'HGNC': '3'})
    st1 = Phosphorylation(a, c, 'T', '185')
    st2 = Dephosphorylation(a, c, 'T', '185')
    st3 = Phosphorylation(None, c, 'T', '185')
    st4 = Dephosphorylation(None, c, 'T', '185')

    return [st1, st2, st3, st4]


@with_model
def test_distinguish_path_polarity1():
    """Test the ability to distinguish a positive from a negative regulation."""
    Monomer('A')
    Monomer('B', ['act'], {'act' : ['y', 'n']})
    Monomer('C', ['T185'], {'T185': ['u', 'p']})
    Parameter('k', 1)
    Rule('A_activate_B', A() + B(act='n') >> A() + B(act='y'), k)
    Rule('B_dephos_C', B(act='y') + C(T185='p') >>
                       B(act='y') + C(T185='u'), k)
    Initial(A(), k)
    Initial(B(act='y'), k)
    Initial(C(T185='p'), k)
    Annotation(A, 'https://identifiers.org/hgnc:1')
    Annotation(B, 'https://identifiers.org/hgnc:2')
    Annotation(C, 'https://identifiers.org/hgnc:3')
    Annotation('A_activate_B', 'A', 'rule_has_subject')
    Annotation('A_activate_B', 'B', 'rule_has_object')
    Annotation('B_dephos_C', 'B', 'rule_has_subject')
    Annotation('B_dephos_C', 'C', 'rule_has_object')
    C.site_annotations = [
            Annotation(('T185', 'p'), 'phosphorylation', 'is_modification'),
            Annotation('T185', 'T', 'is_residue'),
            Annotation('T185', '185', 'is_position'),
        ]
    # Create the model checker
    stmts = _path_polarity_stmt_list()
    mc = PysbModelChecker(model, stmts)
    results = mc.check_model()
    assert len(results) == len(stmts)
    assert isinstance(results[0], tuple)
    path_results = [res[1] for res in results]
    assert path_results[0].paths == []
    assert path_results[1].paths == [(('A_activate_B', 0), ('B_dephos_C', 0),
                                      ('C_T185_p_obs', 1))]
    assert path_results[2].paths == [], path_results[2].paths
    assert path_results[3].paths == \
        [(('B_dephos_C', 0), ('C_T185_p_obs', 1))], path_results[3].paths

@with_model
def test_distinguish_path_polarity2():
    """Test the ability to distinguish a positive from a negative regulation."""
    Monomer('A')
    Monomer('B', ['act'], {'act' : ['y', 'n']})
    Monomer('C', ['T185'], {'T185': ['u', 'p']})
    Parameter('k', 1)
    Rule('A_inhibit_B', A() + B(act='y') >> A() + B(act='n'), k)
    Rule('B_dephos_C', B(act='y') + C(T185='p') >>
                       B(act='y') + C(T185='u'), k)
    Initial(A(), k)
    Initial(B(act='y'), k)
    Initial(C(T185='p'), k)
    Annotation(A, 'https://identifiers.org/hgnc:1')
    Annotation(B, 'https://identifiers.org/hgnc:2')
    Annotation(C, 'https://identifiers.org/hgnc:3')
    Annotation('A_inhibit_B', 'A', 'rule_has_subject')
    Annotation('A_inhibit_B', 'B', 'rule_has_object')
    Annotation('B_dephos_C', 'B', 'rule_has_subject')
    Annotation('B_dephos_C', 'C', 'rule_has_object')
    C.site_annotations = [
            Annotation(('T185', 'p'), 'phosphorylation', 'is_modification'),
            Annotation('T185', 'T', 'is_residue'),
            Annotation('T185', '185', 'is_position'),
        ]
    # Create the model checker
    stmts = _path_polarity_stmt_list()
    mc = PysbModelChecker(model, stmts)
    results = mc.check_model()
    assert len(results) == len(stmts)
    assert isinstance(results[0], tuple)
    assert results[0][1].paths == [(('A_inhibit_B', 0), ('B_dephos_C', 1),
                                    ('C_T185_p_obs', 0))]
    assert results[1][1].paths == []
    assert results[2][1].paths == [(('A_inhibit_B', 0), ('B_dephos_C', 1),
                                    ('C_T185_p_obs', 0))], results[2][1].paths
    assert results[3][1].paths == [(('B_dephos_C', 0), ('C_T185_p_obs', 1))]


def test_check_activation():
    a = Agent('A', db_refs={'HGNC': '1'})
    b = Agent('B', db_refs={'HGNC': '2'})
    c = Agent('C', db_refs={'HGNC': '3'})
    st1 = Activation(a, b)
    st2 = Inhibition(b, c, 'kinase')
    stmts = [st1, st2]
    # Create the model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model(policies='one_step')
    mc = PysbModelChecker(pa.model, stmts)
    results = mc.check_model()
    assert len(results) == len(stmts)
    assert isinstance(results[0], tuple)
    assert results[0][1].paths == [(('A_activates_B_activity', 0),
                                    ('B_activity_active_obs', 0))]
    assert results[1][1].paths == [(('B_deactivates_C_kinase', 0),
                                    ('C_kinase_active_obs', 1))]


def test_check_increase_grounded():
    mmp9 = Agent('MMP9', activity=ActivityCondition('catalytic', True),
                 db_refs={'HGNC': '7176', 'UP': 'P14780'})
    tgfb1 = Agent('TGFB1', db_refs={'HGNC': '11766', 'UP': 'P01137'})
    mapk1 = Agent('MAPK1', db_refs={'HGNC': '6871', 'UP': 'P28482'})
    stmt1 = IncreaseAmount(mmp9, tgfb1)
    stmt2 = IncreaseAmount(mapk1, tgfb1)
    # Make the model out of statement 2 (mapk1), but test with stmt 1 (mmp9)
    pa = PysbAssembler()
    pa.add_statements([stmt2])
    pa.make_model(policies='one_step')
    mc = PysbModelChecker(pa.model, [stmt1])
    results = mc.check_model()
    assert len(results) == 1
    path_result = results[0][1]
    assert path_result.path_found is False
    assert path_result.result_code == 'SUBJECT_MONOMERS_NOT_FOUND'


def test_check_increase_grounded_with_state():
    mmp9_act = Agent('MMP9', activity=ActivityCondition('catalytic', True),
                 db_refs={'HGNC': '7176', 'UP': 'P14780'})
    mmp9 = Agent('MMP9', db_refs={'HGNC': '7176', 'UP': 'P14780'})
    tgfb1 = Agent('TGFB1', db_refs={'HGNC': '11766', 'UP': 'P01137'})
    stmt1 = IncreaseAmount(mmp9_act, tgfb1)
    stmt2 = IncreaseAmount(mmp9, tgfb1)
    # Make the model out of statement 2 (no activity), but test with stmt 1
    # (has activity)
    pa = PysbAssembler()
    pa.add_statements([stmt2])
    pa.make_model(policies='one_step')
    mc = PysbModelChecker(pa.model, [stmt1])
    results = mc.check_model()
    assert len(results) == 1
    path_result = results[0][1]
    assert path_result.path_found is True
    assert stmts_from_pysb_path(
        path_result.paths[0], pa.model, [stmt2]) == [stmt2]


def test_check_activation_grounded():
    mmp9 = Agent('MMP9', activity=ActivityCondition('catalytic', True),
                 db_refs={'HGNC': '7176', 'UP': 'P14780'})
    tgfb1 = Agent('TGFB1', db_refs={'HGNC': '11766', 'UP': 'P01137'})
    mapk1 = Agent('MAPK1', db_refs={'HGNC': '6871', 'UP': 'P28482'})
    stmt1 = Activation(mmp9, tgfb1)
    stmt2 = Activation(mapk1, tgfb1)
    # Make the model out of statement 2 (mapk1), but test with stmt 1 (mmp9)
    pa = PysbAssembler()
    pa.add_statements([stmt2])
    pa.make_model(policies='one_step')
    mc = PysbModelChecker(pa.model, [stmt1])
    results = mc.check_model()
    assert len(results) == 1
    path_result = results[0][1]
    assert path_result.path_found is False
    assert path_result.result_code == 'SUBJECT_MONOMERS_NOT_FOUND'


@with_model
def test_none_phosphorylation_stmt():
    # Create the statement
    b = Agent('B', db_refs={'HGNC': '2'})
    st1 = Phosphorylation(None, b, 'T', '185')
    st2 = Phosphorylation(None, b, 'Y', '187')
    stmts = [st1, st2]
    # Now create the PySB model
    Monomer('A')
    Monomer('B', ['T185', 'Y187'], {'T185':['u', 'p'], 'Y187': ['u', 'p']})
    Rule('A_phos_B', A() + B(T185='u') >> A() + B(T185='p'),
         Parameter('k', 1))
    Initial(A(), Parameter('A_0', 100))
    Initial(B(T185='u', Y187='p'), Parameter('B_0', 100))
    Annotation(A, 'https://identifiers.org/hgnc:1')
    Annotation(B, 'https://identifiers.org/hgnc:2')
    B.site_annotations = [
        Annotation(('T185', 'p'), 'phosphorylation', 'is_modification'),
        Annotation('T185', 'T', 'is_residue'),
        Annotation('T185', '185', 'is_position'),
        Annotation(('Y187', 'p'), 'phosphorylation', 'is_modification'),
        Annotation('Y187', 'Y', 'is_residue'),
        Annotation('Y187', '187', 'is_position'),
    ]
    mc = PysbModelChecker(model, stmts)
    results = mc.check_model()
    assert len(results) == 2
    assert isinstance(results[0], tuple)
    assert results[0][0] == st1
    assert results[0][1].paths == [(('A_phos_B', 0), ('B_T185_p_obs', 0))]
    assert results[1][0] == st2
    assert results[1][1].paths == []


@with_model
def test_phosphorylation_annotations():
    # Create the statement
    a = Agent('MEK1', db_refs={'HGNC': '6840'})
    b = Agent('ERK2', db_refs={'HGNC': '6871'})
    st1 = Phosphorylation(a, b, 'T', '185')
    st2 = Phosphorylation(a, b, None, None)
    st3 = Phosphorylation(a, b, 'Y', '187')
    # Now create the PySB model
    Monomer('A_monomer')
    Monomer('B_monomer', ['Thr185', 'Y187'],
            {'Thr185': ['un', 'phos'], 'Y187': ['u', 'p']})
    Rule('A_phos_B', A_monomer() + B_monomer(Thr185='un') >>
                     A_monomer() + B_monomer(Thr185='phos'),
         Parameter('k', 1))
    Initial(A_monomer(), Parameter('A_0', 100))
    Initial(B_monomer(Thr185='un', Y187='u'), Parameter('B_0', 100))
    # Add agent grounding
    Annotation(A_monomer, 'https://identifiers.org/hgnc:6840')
    Annotation(B_monomer, 'https://identifiers.org/hgnc:6871')
    Annotation('A_phos_B', 'A_monomer', 'rule_has_subject')
    Annotation('A_phos_B', 'B_monomer', 'rule_has_object')
    # Add annotations to the sites/states of the Monomer itself
    B_annot = [
        Annotation('Thr185', 'T', 'is_residue'),
        Annotation('Thr185', '185', 'is_position'),
        Annotation(('Thr185', 'phos'), 'phosphorylation', 'is_modification'),
        Annotation('Y187', 'Y', 'is_residue'),
        Annotation('Y187', '187', 'is_position'),
        Annotation(('Y187', 'p'), 'phosphorylation', 'is_modification'),
    ]
    B_monomer.site_annotations = B_annot
    mc = PysbModelChecker(model, [st1, st2, st3])
    results = mc.check_model()
    assert len(results) == 3
    assert isinstance(results[0], tuple)
    assert results[0][0] == st1
    assert results[0][1].paths == [(('A_phos_B', 0),
                                    ('B_monomer_Thr185_phos_obs', 0))]
    assert results[1][0] == st2
    assert results[1][1].paths == [(('A_phos_B', 0),
                                    ('B_monomer_Thr185_phos_obs', 0))]
    assert results[2][0] == st3
    assert results[2][1].paths == []


@with_model
def test_activation_annotations():
    # Create the statement
    a = Agent('MEK1', db_refs={'HGNC': '6840'})
    b = Agent('ERK2', db_refs={'HGNC': '6871'})
    st1 = Phosphorylation(a, b, 'T', '185')
    st2 = Activation(a, b)
    st3 = Activation(a, b, 'kinase')
    # Now create the PySB model
    Monomer('A_monomer')
    Monomer('B_monomer', ['Thr185', 'Y187'],
            {'Thr185': ['un', 'phos'], 'Y187': ['u', 'p']})
    Rule('A_phos_B', A_monomer() + B_monomer(Thr185='un') >>
                     A_monomer() + B_monomer(Thr185='phos'),
         Parameter('k', 1))
    Initial(A_monomer(), Parameter('A_0', 100))
    Initial(B_monomer(Thr185='un', Y187='u'), Parameter('B_0', 100))
    # Add agent grounding
    Annotation(A_monomer, 'https://identifiers.org/hgnc:6840')
    Annotation(B_monomer, 'https://identifiers.org/hgnc:6871')
    Annotation(B_monomer, {'Thr185':'phos'}, 'has_active_pattern')
    Annotation('A_phos_B', 'A_monomer', 'rule_has_subject')
    Annotation('A_phos_B', 'B_monomer', 'rule_has_object')
    # Add annotations to the sites/states of the Monomer itself
    B_annot = [
        Annotation('Thr185', 'T', 'is_residue'),
        Annotation('Thr185', '185', 'is_position'),
        Annotation(('Thr185', 'phos'), 'phosphorylation', 'is_modification'),
        Annotation('Y187', 'Y', 'is_residue'),
        Annotation('Y187', '187', 'is_position'),
        Annotation(('Y187', 'p'), 'phosphorylation', 'is_modification'),
    ]
    B_monomer.site_annotations = B_annot
    mc = PysbModelChecker(model, [st1, st2, st3])
    results = mc.check_model()
    assert len(results) == 3
    assert isinstance(results[0], tuple)
    assert results[0][0] == st1
    assert results[0][1].paths == [(('A_phos_B', 0),
                                    ('B_monomer_Thr185_phos_obs', 0))]
    assert results[1][0] == st2
    assert results[1][1].paths == [(('A_phos_B', 0),
                                    ('B_monomer_Thr185_phos_obs', 0))]
    assert results[2][0] == st3
    assert results[1][1].paths == [(('A_phos_B', 0),
                                    ('B_monomer_Thr185_phos_obs', 0))]



def test_multitype_path():
    """Test causal chain involving Complex, Gef, Activation"""
    egfr = Agent('EGFR', db_refs={'HGNC': '3236'})
    grb2 = Agent('GRB2', db_refs={'HGNC': '4566'})
    grb2_egfr = Agent('GRB2', bound_conditions=[BoundCondition(egfr)],
                      db_refs={'HGNC': '4566'})
    sos1 = Agent('SOS1', db_refs={'HGNC': '11187'}, )
    sos1_grb2 = Agent('SOS1', bound_conditions=[BoundCondition(grb2)],
                      db_refs={'HGNC': '11187'}, )
    kras = Agent('KRAS', db_refs={'HGNC': '6407'})
    kras_g = Agent('KRAS', activity=ActivityCondition('gtpbound', True),
                   db_refs={'HGNC': '6407'})
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})

    def check_stmts(stmts, paths):
        pa = PysbAssembler()
        pa.add_statements(stmts)
        pa.make_model(policies='one_step')
        stmts_to_check = [
                Activation(egfr, kras, 'gtpbound'),
                Activation(egfr, braf, 'kinase')
            ]
        mc = PysbModelChecker(pa.model, stmts_to_check)
        results = mc.check_model()
        assert len(results) == len(stmts_to_check)
        assert isinstance(results[0], tuple)
        assert results[0][1].paths == paths[0], results[0][1].paths
        assert results[1][1].paths == paths[1], results[1][1].paths
    # Check with the ActiveForm
    stmts1 = [
        Complex([egfr, grb2]),
        Complex([sos1, grb2_egfr]),
        ActiveForm(sos1_grb2, 'activity', True),
        Activation(sos1_grb2, kras, 'gtpbound'),
        Activation(kras_g, braf, 'kinase')
      ]
    check_stmts(stmts1, ([(('EGFR_GRB2_bind', 0), ('SOS1_GRB2_EGFR_bind', 0),
                           ('SOS1_GRB2_activates_KRAS_gtpbound', 0),
                           ('KRAS_gtpbound_active_obs', 0))],
                         [(('EGFR_GRB2_bind', 0), ('SOS1_GRB2_EGFR_bind', 0),
                           ('SOS1_GRB2_activates_KRAS_gtpbound', 0),
                           ('KRAS_gtp_activates_BRAF_kinase', 0),
                           ('BRAF_kinase_active_obs', 0))]))
    # Check without the ActiveForm
    stmts2 = [
        Complex([egfr, grb2]),
        Complex([sos1, grb2_egfr]),
        Gef(sos1_grb2, kras),
        Activation(kras_g, braf, 'kinase')
      ]
    check_stmts(stmts2, ([(('EGFR_GRB2_bind', 0), ('SOS1_GRB2_EGFR_bind', 0),
                           ('SOS1_GRB2_activates_KRAS', 0),
                           ('KRAS_gtpbound_active_obs', 0))],
                         [(('EGFR_GRB2_bind', 0), ('SOS1_GRB2_EGFR_bind', 0),
                           ('SOS1_GRB2_activates_KRAS', 0),
                           ('KRAS_gtp_activates_BRAF_kinase', 0),
                           ('BRAF_kinase_active_obs', 0))]))


def test_grounded_modified_enzyme():
    """Check if the model checker can use semantic annotations to match mods
    on the enzyme, not just the substrate, of a phosphorylation statement."""
    mek_s202 = Agent('MEK1', mods=[ModCondition('phosphorylation', 'S', '202')],
                     db_refs={'HGNC': '6840'})
    mek_phos = Agent('MEK1', mods=[ModCondition('phosphorylation', None, None)],
                     db_refs={'HGNC': '6840'})
    erk = Agent('ERK2', db_refs={'HGNC': '6871'})
    stmt_to_model = Phosphorylation(mek_s202, erk, None, None)
    stmt_to_check = Phosphorylation(mek_phos, erk, None, None)
    pa = PysbAssembler()
    pa.add_statements([stmt_to_model])
    pa.make_model(policies='one_step')
    mc = PysbModelChecker(pa.model, [stmt_to_check])
    results = mc.check_model()
    assert len(results) == 1
    assert results[0][0] == stmt_to_check
    assert results[0][1].paths == \
        [(('MEK1_phosphoS202_phosphorylation_ERK2_phospho', 0),
          ('ERK2_phospho_p_obs', 0))]
    assert stmts_from_pysb_path(
        results[0][1].paths[0], pa.model, [stmt_to_model]) == [stmt_to_model]


def test_check_ubiquitination():
    xiap = Agent('XIAP', db_refs={'HGNC': '592'})
    casp3 = Agent('CASP3', db_refs={'HGNC': '1504'})
    stmt = Ubiquitination(xiap, casp3)
    pysba = PysbAssembler()
    pysba.add_statements([stmt])
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [stmt])
    checks = mc.check_model()
    assert len(checks) == 1
    assert isinstance(checks[0], tuple)
    assert checks[0][0] == stmt
    assert checks[0][1].paths == [(('XIAP_ubiquitination_CASP3_ub', 0),
                                   ('CASP3_ub_y_obs', 0))]
    assert stmts_from_pysb_path(
        checks[0][1].paths[0], pysba.model, [stmt]) == [stmt]


def test_check_rule_subject1():
    mek = Agent('MEK1', db_refs={'HGNC': '6840'})
    erk = Agent('ERK2', db_refs={'HGNC': '6871'})
    stmt = Phosphorylation(mek, erk)
    pysba = PysbAssembler()
    pysba.add_statements([stmt])
    pysba.make_model(policies='one_step')
    # Check against stmt: should not validate ERK phosphorylates ERK
    stmt_to_check = Phosphorylation(erk, erk)
    mc = PysbModelChecker(pysba.model, [stmt_to_check])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == stmt_to_check
    assert checks[0][1].paths == [], checks[0][1].paths


def test_gef_activation():
    sos = Agent('SOS1', db_refs={'HGNC': '1'})
    ras = Agent('KRAS', db_refs={'HGNC': '2'})
    gef_stmt = Gef(sos, ras)
    act_stmt = Activation(sos, ras, 'gtpbound')
    # Check that the activation is satisfied by the Gef
    pysba = PysbAssembler()
    pysba.add_statements([gef_stmt])
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [act_stmt])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == act_stmt
    assert checks[0][1].paths == [(('SOS1_activates_KRAS', 0),
                                   ('KRAS_gtpbound_active_obs', 0))]
    assert stmts_from_pysb_path(
        checks[0][1].paths[0], pysba.model, [gef_stmt]) == [gef_stmt]
    # TODO TODO TODO
    """
    # Check that the Gef is satisfied by the Activation
    # This currently doesn't work because Gef statements aren't checked
    pysba = PysbAssembler()
    pysba.add_statements([act_stmt])
    pysba.make_model(policies='one_step')
    mc = ModelChecker(pysba.model, [gef_stmt])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == gef_stmt
    assert checks[0][1] == True
    """


def test_gef_rasgtp():
    sos = Agent('SOS1', db_refs={'HGNC': '1'})
    ras = Agent('KRAS', db_refs={'HGNC': '2'})
    ras_gtp = Agent('KRAS', activity=ActivityCondition('gtpbound', True),
                    db_refs={'HGNC': '2'})
    raf = Agent('BRAF', db_refs={'HGNC': '3'})
    gef_stmt = Gef(sos, ras)
    rasgtp_stmt = GtpActivation(ras_gtp, raf, 'kinase')
    act_stmt = Activation(sos, raf, 'kinase')
    # Check that the activation is satisfied by the Gef
    pysba = PysbAssembler()
    pysba.add_statements([gef_stmt, rasgtp_stmt])
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [act_stmt])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == act_stmt
    assert checks[0][1].paths == [(('SOS1_activates_KRAS', 0),
                                   ('KRAS_gtp_activates_BRAF_kinase', 0),
                                   ('BRAF_kinase_active_obs', 0))], \
        checks[0][1].paths
    assert stmts_from_pysb_path(checks[0][1].paths[0], pysba.model, [
        gef_stmt, rasgtp_stmt]) == [gef_stmt, rasgtp_stmt]


def test_gef_rasgtp_phos():
    sos = Agent('SOS1', db_refs={'HGNC': '1'})
    ras = Agent('KRAS', db_refs={'HGNC': '2'})
    ras_a = Agent('KRAS', activity=ActivityCondition('gtpbound', True),
                  db_refs={'HGNC': '2'})
    raf = Agent('BRAF', db_refs={'HGNC': '3'})
    raf_a = Agent('BRAF', activity=ActivityCondition('kinase', True),
                  db_refs={'HGNC': '3'})
    mek = Agent('MEK', db_refs={'HGNC': '4'})
    gef_stmt = Gef(sos, ras)
    rasgtp_stmt = GtpActivation(ras_a, raf, 'kinase')
    phos = Phosphorylation(raf_a, mek)
    stmt_to_check = Phosphorylation(sos, mek)
    # Assemble and check
    pysba = PysbAssembler()
    pysba.add_statements([gef_stmt, rasgtp_stmt, phos])
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [stmt_to_check])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == stmt_to_check
    assert checks[0][1].paths == [(('SOS1_activates_KRAS', 0),
                                   ('KRAS_gtp_activates_BRAF_kinase', 0),
                                   ('BRAF_kin_phosphorylation_MEK_phospho', 0),
                                   ('MEK_phospho_p_obs', 0))], \
        checks[0][1].paths
    assert stmts_from_pysb_path(checks[0][1].paths[0], pysba.model, [
        gef_stmt, rasgtp_stmt, phos]) == [gef_stmt, rasgtp_stmt, phos]


def test_gap_activation():
    nf1 = Agent('NF1', db_refs={'HGNC': '1'})
    ras = Agent('KRAS', db_refs={'HGNC': '2'})
    gap_stmt = Gap(nf1, ras)
    act_stmt = Inhibition(nf1, ras, 'gtpbound')
    # Check that the activation is satisfied by the Gap
    pysba = PysbAssembler()
    pysba.add_statements([gap_stmt])
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [act_stmt])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == act_stmt
    assert checks[0][1].paths == [(('NF1_deactivates_KRAS', 0),
                                   ('KRAS_gtpbound_active_obs', 1))]
    assert stmts_from_pysb_path(
        checks[0][1].paths[0], pysba.model, [gap_stmt]) == [gap_stmt]
    # TODO TODO TODO
    """
    # Check that the Gap is satisfied by the Activation
    # This currently doesn't work because Gap statements aren't checked by
    # the ModelChecker
    pysba = PysbAssembler()
    pysba.add_statements([act_stmt])
    pysba.make_model(policies='one_step')
    mc = ModelChecker(pysba.model, [gap_stmt])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == gap_stmt
    assert checks[0][1] == True
    """


def test_gap_rasgtp():
    nf1 = Agent('NF1', db_refs={'HGNC': '1'})
    ras = Agent('KRAS', db_refs={'HGNC': '2'})
    ras_g = Agent('KRAS', activity=ActivityCondition('gtpbound', True),
                  db_refs={'HGNC': '2'})
    raf = Agent('BRAF', db_refs={'HGNC': '3'})
    gap_stmt = Gap(nf1, ras)
    rasgtp_stmt = GtpActivation(ras_g, raf, 'kinase')
    act_stmt = Inhibition(nf1, raf, 'kinase')
    # Check that the activation is satisfied by the Gap
    pysba = PysbAssembler()
    pysba.add_statements([gap_stmt, rasgtp_stmt])
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [act_stmt])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == act_stmt
    assert checks[0][1].paths == [(('NF1_deactivates_KRAS', 0),
                                   ('KRAS_gtp_activates_BRAF_kinase', 1),
                                   ('BRAF_kinase_active_obs', 1))], \
        checks[0][1].paths
    assert stmts_from_pysb_path(checks[0][1].paths[0], pysba.model, [
        gap_stmt, rasgtp_stmt]) == [gap_stmt, rasgtp_stmt]


def test_gap_rasgtp_phos():
    nf1 = Agent('NF1', db_refs={'HGNC': '1'})
    ras = Agent('KRAS', db_refs={'HGNC': '2'})
    ras_g = Agent('KRAS', activity=ActivityCondition('gtpbound', True),
                  db_refs={'HGNC': '2'})
    raf = Agent('BRAF', db_refs={'HGNC': '3'})
    raf_a = Agent('BRAF', activity=ActivityCondition('kinase', True),
                  db_refs={'HGNC': '3'})
    mek = Agent('MEK', db_refs={'HGNC': '4'})
    gap_stmt = Gap(nf1, ras)
    rasgtp_stmt = GtpActivation(ras_g, raf, 'kinase')
    phos = Phosphorylation(raf_a, mek)
    stmt_to_check = Dephosphorylation(nf1, mek)
    # Assemble and check
    pysba = PysbAssembler()
    pysba.add_statements([gap_stmt, rasgtp_stmt, phos])
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [stmt_to_check])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == stmt_to_check
    assert checks[0][1].paths == \
        [(('NF1_deactivates_KRAS', 0),
          ('KRAS_gtp_activates_BRAF_kinase', 1),
          ('BRAF_kin_phosphorylation_MEK_phospho', 1),
          ('MEK_phospho_p_obs', 1))], checks[0][1].paths
    assert stmts_from_pysb_path(checks[0][1].paths[0], pysba.model, [
        gap_stmt, rasgtp_stmt, phos]) == [gap_stmt, rasgtp_stmt, phos]


def test_increase_amount():
    tp53 = Agent('TP53', db_refs={'HGNC': '1'})
    x = Agent('X', db_refs={'HGNC': 2})
    mdm2 = Agent('MDM2', db_refs={'HGNC': '3'})
    stmts = [IncreaseAmount(tp53, x), IncreaseAmount(x, mdm2)]
    stmt_to_check = IncreaseAmount(tp53, mdm2)
    pysba = PysbAssembler()
    pysba.add_statements(stmts)
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [stmt_to_check])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == stmt_to_check
    assert checks[0][1].paths == [(('TP53_produces_X', 0),
                                   ('X_produces_MDM2', 0),
                                   ('MDM2__obs', 0))]
    assert stmts_from_pysb_path(
        checks[0][1].paths[0], pysba.model, stmts) == stmts


def test_decrease_amount():
    tp53 = Agent('TP53', db_refs={'HGNC': '1'})
    tp53u = Agent('TP53', mods=[ModCondition('ubiquitination')],
                  db_refs={'HGNC': '1'})
    mdm2 = Agent('MDM2', db_refs={'HGNC': '3'})
    stmts = [IncreaseAmount(tp53, mdm2),
             Ubiquitination(mdm2, tp53), DecreaseAmount(None, tp53u)]
    stmt_to_check = DecreaseAmount(tp53, tp53)
    pysba = PysbAssembler()
    pysba.add_statements(stmts)
    pysba.make_model(policies='one_step')
    mc = PysbModelChecker(pysba.model, [stmt_to_check])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == stmt_to_check
    assert checks[0][1].paths == [(('TP53_produces_MDM2', 0),
                                   ('MDM2_ubiquitination_TP53_ub', 0),
                                   ('TP53_ub_degraded', 0),
                                   ('TP53__obs', 1))]
    assert stmts_from_pysb_path(
        checks[0][1].paths[0], pysba.model, stmts) == stmts


def test_stmt_from_rule():
    mek = Agent('MEK1', db_refs={'HGNC': '6840'})
    erk = Agent('ERK2', db_refs={'HGNC': '6871'})
    st = Phosphorylation(mek, erk, 'T', '185')
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    rule_name = pa.model.rules[0].name
    stmt = stmt_from_rule(rule_name, pa.model, [st])
    assert stmt == st


def test_activate_via_mod():
    mek = Agent('MEK1', db_refs={'HGNC': '6840'})
    erk = Agent('ERK2', db_refs={'HGNC': '6871'})
    erka = Agent('ERK2', mods=[ModCondition('phosphorylation', 'T', '185')],
                 db_refs={'HGNC': '6871'})
    st1 = Phosphorylation(mek, erk, 'T', '185')
    st2 = ActiveForm(erka, 'activity', True)
    st3 = Activation(mek, erk)
    pa = PysbAssembler()
    pa.add_statements([st1, st2])
    pa.make_model()
    mc = PysbModelChecker(pa.model, [st3])
    checks = mc.check_model()
    # Make sure it checks out to True
    assert checks[0][1].path_found


def test_observables():
    mek = Agent('MEK1', db_refs={'HGNC': '6840'})
    erk = Agent('ERK2', db_refs={'HGNC': '6871'})
    erkp = Agent('ERK2', mods=[ModCondition('phosphorylation', 'T', '185')],
                 db_refs={'HGNC': '6871'})
    st1 = Phosphorylation(mek, erk, 'T', '185')
    st2 = ActiveForm(erkp, 'activity', True)
    st3 = Activation(mek, erk)
    pa = PysbAssembler()
    pa.add_statements([st1, st2])
    pa.make_model()
    mc = PysbModelChecker(pa.model, [st1, st3], agent_obs=[erkp])
    checks = mc.check_model()
    assert checks[0][1].path_found
    assert checks[1][1].path_found
    # Only 1 observable should be created
    assert len(mc.model.observables) == 1


"""
def test_check_rule_subject_bound_condition():
    braf = Agent('BRAF', db_refs={'HGNC': '1'})
    raf1 = Agent('RAF1', db_refs={'HGNC': '2'})
    braf_raf1 = Agent('BRAF', bound_conditions=[BoundCondition(raf1)],
                      db_refs={'HGNC': '1'})
    mek = Agent('MEK1', db_refs={'HGNC': '6840'})
    stmt = Phosphorylation(braf_raf1, mek)
    pysba = PysbAssembler()
    pysba.add_statements([stmt])
    pysba.make_model(policies='one_step')
    # Check against stmt: should indicate that RAF1 is causally linked to MEK
    # phosphorylation
    stmt_to_check = Phosphorylation(raf1, mek)
    mc = ModelChecker(pysba.model, [stmt_to_check])
    checks = mc.check_model()
    assert len(checks) == 1
    assert checks[0][0] == stmt_to_check
    assert checks[0][1] == True

def test_activation_subtype():
    sos1 = Agent('SOS1', db_refs={'HGNC':'11187'})
    kras = Agent('KRAS', db_refs={'HGNC':'6407'})
    stmts = [Activation(sos1, kras, 'gtpbound')]
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model(policies='one_step')
    stmts_to_check = [Activation(sos1, kras, 'activity')]
    mc = ModelChecker(pa.model, stmts_to_check)
    results = mc.check_model()
    assert len(results) == len(stmts_to_check)
    assert isinstance(results[0], tuple)
    assert results[0][1] == True

def test_check_autophosphorylation():
    egfr = Agent('EGFR', db_refs={'HGNC':'3236'})
    stmts = [Autophosphorylation(egfr, 'Y', '1016')]
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model(policies='one_step')
    stmts_to_check = [Phosphorylation(None, egfr),
                      Phosphorylation(None, egfr, 'Y', '1016')]
    mc = ModelChecker(pa.model, stmts_to_check)
    results = mc.check_model()
    assert len(results) == len(stmts_to_check)
    assert isinstance(results[0], tuple)
    assert results[0][1] == True
    assert results[1][1] == True

def test_check_transphosphorylation():
    egfr = Agent('EGFR', db_refs={'HGNC':'3236'})
    erbb2_egfr = Agent('ERBB2', bound_conditions=[BoundCondition(egfr)],
                       db_refs={'HGNC':'3430'})
    stmts = [Transphosphorylation(erbb2_egfr, 'Y', '1016')]
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model(policies='one_step')
    stmts_to_check = [Phosphorylation(None, egfr),
                      Phosphorylation(None, egfr, 'Y', '1016')]
    mc = ModelChecker(pa.model, stmts_to_check)
    results = mc.check_model()
    assert len(results) == len(stmts_to_check)
    assert isinstance(results[0], tuple)
    assert results[0][1] == True
    assert results[1][1] == True

"""


@unittest.skip('Skip path score tests for now')
def test_model_check_data():
    # Create a set of statements
    a = Agent('A', db_refs={'HGNC': '1'})
    b = Agent('B', db_refs={'HGNC': '2'})
    b_phos = Agent('B', mods=[ModCondition('phosphorylation')],
                   db_refs={'HGNC': '2'})
    c = Agent('C', db_refs={'HGNC': '3'})
    c_phos = Agent('C', mods=[ModCondition('phosphorylation')],
                   db_refs={'HGNC': '3'})
    d = Agent('D', db_refs={'HGNC': '4'})
    d_phos = Agent('D', mods=[ModCondition('phosphorylation')],
                   db_refs={'HGNC': '4'})

    # Two paths from A to D: One going through B and another through C
    st1 = Phosphorylation(a, b)
    st2 = Phosphorylation(b_phos, d)
    st3 = Phosphorylation(a, c)
    st4 = Phosphorylation(c_phos, d)
    # Statements/Data agents for checking
    stmt_to_check = Phosphorylation(a, d)
    agent_obs = [b_phos, c_phos, d_phos]
    # Make model
    pa = PysbAssembler()
    pa.add_statements([st1, st2, st3, st4])
    pa.make_model(policies='one_step')
    mc = PysbModelChecker(pa.model, [stmt_to_check], agent_obs)
    results = mc.check_model(max_paths=5)
    # Create observable
    assert len(results) == 1
    pr = results[0][1]
    res = pr.paths[0:2]
    assert len(res) == 2
    p1 = (('A_phosphorylation_B_phospho', 0),
          ('B_phospho_phosphorylation_D_phospho', 0),
          ('D_phospho_p_obs', 0))
    assert p1 in res
    p2 = (('A_phosphorylation_C_phospho', 0),
          ('C_phospho_phosphorylation_D_phospho', 0),
          ('D_phospho_p_obs', 0))
    assert p2 in res
    # Now, a vector linking agents with values, expressed at first as
    # +/- 1
    # This data should ensure that the path through B should be more highly
    # ranked than the path through C
    data = {b_phos: 1, c_phos: -1, d_phos: 1}
    paths = results[0][1].paths
    scored_paths = mc.score_paths(paths, data)
    assert scored_paths[0][0] == p1
    assert scored_paths[1][0] == p2
    assert scored_paths[0][1] > scored_paths[1][1]


def test_prune_influence_map():
    kin = Agent('Kinase', db_refs={'HGNC': '1'})
    phos = Agent('Phosphatase', db_refs={'HGNC': '2'})
    subs = Agent('Substrate', db_refs={'HGNC': '3'})
    st1 = Phosphorylation(kin, subs)
    st2 = Dephosphorylation(phos, subs)
    pa = PysbAssembler()
    pa.add_statements([st1, st2])
    pa.make_model(policies='one_step')
    mc = PysbModelChecker(pa.model, [st1])
    im = mc.get_im()
    remove_im_params(pa.model, im)
    mc.prune_influence_map()
    im = mc.get_im()
    assert len(im.nodes()) == 3
    assert len(im.edges()) == 2
    # Smoke test to make sure drawing works
    mc.draw_im(os.devnull)


def test_prune_influence_map_subj_obj():
    def ag(gene_name, delta=None):
        return Event(Agent(gene_name,
                           db_refs={'HGNC':
                                    hgnc_client.get_hgnc_id(gene_name)}),
                     delta=delta)
    mek = ag('MAP2K1')
    erk = ag('MAPK1')
    erk_neg = ag('MAPK1', QualitativeDelta(polarity=-1))
    mek2 = ag('MAP2K2')
    mek2_neg = ag('MAP2K2', QualitativeDelta(polarity=-1))

    s1 = Influence(mek, erk)
    s2 = Influence(mek2, erk_neg)
    s3 = Influence(erk, mek2_neg)
    # To check:
    s4 = Influence(mek, mek2)
    # Make the model
    pa = PysbAssembler()
    pa.add_statements([s1, s2, s3])
    model = pa.make_model()
    # Check the model
    mc = PysbModelChecker(model, [s4])
    pr_before = mc.check_statement(s4)
    assert pr_before.result_code == 'PATHS_FOUND', pr_before
    # Now prune the influence map
    mc.graph = None
    mc.get_graph(prune_im=True, prune_im_degrade=True, prune_im_subj_obj=True)
    pr_after = mc.check_statement(s4)
    assert pr_after.result_code == 'NO_PATHS_FOUND', pr_after


def test_prune_influence_map_degrade_bind():
    deg = DecreaseAmount(None, Agent('X'))
    bind = Complex([Agent('X'), Agent('Y')])
    pa = PysbAssembler([deg, bind])
    model = pa.make_model()
    # Check the model
    mc = PysbModelChecker(model)
    mc.prune_influence_map()
    im = mc.get_im()
    assert len(im.edges()) == 3, im.edges()
    mc.prune_influence_map_degrade_bind_positive([deg, bind])
    im = mc.get_im()
    assert len(im.edges()) == 2, im.edges()


def test_add_namespaces():
    a = Agent('a', db_refs={'HGNC': '1'})
    b = Agent('b', db_refs={'CHEBI': '2'})
    c = Agent('c', db_refs={'GO': '3'})
    d = Agent('d', db_refs={'HGNC': '4'})
    stmts = [Activation(a, b), Inhibition(b, c), IncreaseAmount(a, c),
             DecreaseAmount(c, d)]
    pa = PysbAssembler(stmts)
    model = pa.make_model()
    mc = PysbModelChecker(
        model, statements=[IncreaseAmount(a, d)], model_stmts=stmts)
    mc.get_graph(add_namespaces=True)
    im = mc.get_im()
    # All nodes have ns
    assert {'a_activates_b_activity', 'a_produces_c', 'd__obs',
            'b_deactivates_c_activity', 'c_degrades_d'} == set(im.nodes)
    for n, data in im.nodes(data=True):
        assert data.get('ns')
        if n in ['a_activates_b_activity', 'a_produces_c', 'd__obs']:
            assert data['ns'] =='HGNC'
        elif n == 'b_deactivates_c_activity':
            assert data['ns'] == 'CHEBI'
        elif n == 'c_degrades_d':
            assert data['ns'] == 'GO'


@unittest.skip('Skip sampling tests for now')
def test_weighted_sampling1():
    """Test sampling with different path lengths but no data."""
    os.environ['TEST_FLAG'] = 'TRUE'
    mc = ModCondition('phosphorylation')
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    map2k1_phos = Agent('MAP2K1', mods=[mc], db_refs={'HGNC': '6840'})
    mapk1 = Agent('MAPK1', db_refs={'HGNC': '6871'})
    mapk1_phos = Agent('MAPK1', mods=[mc], db_refs={'HGNC': '6871'})
    jun = Agent('JUN', db_refs={'HGNC': '6204'})
    stmt_to_check = Phosphorylation(braf, jun)
    stmts = [stmt_to_check,
             Phosphorylation(braf, map2k1),
             Phosphorylation(map2k1_phos, jun),
             Phosphorylation(map2k1_phos, mapk1),
             Phosphorylation(mapk1_phos, jun)]
    # Make model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model(policies='one_step')
    # Make the model checker and prune the influence map
    mc = PysbModelChecker(pa.model, [stmt_to_check], do_sampling=True, seed=1)
    mc.prune_influence_map()
    # Seed the random number generator
    np.random.seed(1)
    results = mc.check_model(max_path_length=5, max_paths=100)
    assert type(results) == list
    assert len(results) == 1
    stmt_tuple = results[0]
    assert len(stmt_tuple) == 2
    assert stmt_tuple[0] == stmt_to_check
    path_result = stmt_tuple[1]
    assert type(path_result) == PathResult
    path_lengths = [len(p) for p in path_result.paths]
    assert max(path_lengths) <= 5
    # There are two distinct paths
    assert len(set(path_result.paths)) == 3
    path_ctr = Counter(path_result.paths)
    assert path_ctr[(('BRAF_phosphorylation_JUN_phospho', 1),
                     ('JUN_phospho_p_obs', 1))] == 46, \
        path_ctr[(('BRAF_phosphorylation_JUN_phospho', 1),
                  ('JUN_phospho_p_obs', 1))]
    assert path_ctr[(('BRAF_phosphorylation_MAP2K1_phospho', 1),
                     ('MAP2K1_phospho_phosphorylation_JUN_phospho', 1),
                     ('JUN_phospho_p_obs', 1))] == 22, path_ctr
    assert path_ctr[(('BRAF_phosphorylation_MAP2K1_phospho', 1),
                     ('MAP2K1_phospho_phosphorylation_MAPK1_phospho', 1),
                     ('MAPK1_phospho_phosphorylation_JUN_phospho', 1),
                     ('JUN_phospho_p_obs', 1))] == 32, path_ctr


@unittest.skip('Skip sampling tests for now')
def test_weighted_sampling2():
    """Test sampling with abundances but no tail probabilities from data."""
    os.environ['TEST_FLAG'] = 'TRUE'
    map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    mapk1 = Agent('MAPK1', db_refs={'HGNC': '6871'})
    mapk3 = Agent('MAPK3', db_refs={'HGNC': '6877'})
    mc = ModCondition('phosphorylation')
    mapk1_phos = Agent('MAPK1', mods=[mc], db_refs={'HGNC': '6871'})
    mapk3_phos = Agent('MAPK3', mods=[mc], db_refs={'HGNC': '6877'})
    jun = Agent('JUN', db_refs={'HGNC': '6204'})
    st1 = Phosphorylation(map2k1, mapk1)
    st2 = Phosphorylation(map2k1, mapk3)
    st3 = Phosphorylation(mapk1_phos, jun)
    st4 = Phosphorylation(mapk3_phos, jun)
    stmt_to_check = Phosphorylation(map2k1, jun)
    # Make model
    pa = PysbAssembler()
    pa.add_statements([st1, st2, st3, st4])
    pa.make_model(policies='one_step')
    # Set the initial conditions
    mapk1_monomer = pa.model.all_components()['MAPK1']
    mapk3_monomer = pa.model.all_components()['MAPK3']
    set_base_initial_condition(pa.model, mapk1_monomer, 75)
    set_base_initial_condition(pa.model, mapk3_monomer, 25)
    # Make the model checker and prune the influence map
    # Setting do_sampling to False should yield the default enumeration
    # behavior
    mc = PysbModelChecker(pa.model, [stmt_to_check], do_sampling=False)
    mc.prune_influence_map()
    results = mc.check_model(max_paths=5)
    path_result = results[0][1]
    assert len(path_result.paths) == 2
    enum_paths = path_result.paths
    # Now, try sampling
    mc = PysbModelChecker(pa.model, [stmt_to_check], do_sampling=True, seed=1)
    mc.prune_influence_map()
    results = mc.check_model(max_path_length=5, max_paths=1000)
    assert type(results) == list
    assert len(results) == 1
    stmt_tuple = results[0]
    assert len(stmt_tuple) == 2
    assert stmt_tuple[0] == stmt_to_check
    path_result = stmt_tuple[1]
    assert type(path_result) == PathResult
    path_lengths = [len(p) for p in path_result.paths]
    assert max(path_lengths) <= 5
    # There are two distinct paths
    assert set(enum_paths) == set(path_result.paths)
    path_ctr = Counter(path_result.paths)
    mapk1_count = path_ctr[(('MAP2K1_phosphorylation_MAPK1_phospho', 1),
                            ('MAPK1_phospho_phosphorylation_JUN_phospho', 1),
                            ('JUN_phospho_p_obs', 1))]
    mapk3_count = path_ctr[(('MAP2K1_phosphorylation_MAPK3_phospho', 1),
                            ('MAPK3_phospho_phosphorylation_JUN_phospho', 1),
                            ('JUN_phospho_p_obs', 1))]
    assert mapk1_count == 750, mapk1_count
    assert mapk3_count == 250, mapk3_count


@unittest.skip('Skip sampling tests for now')
def test_weighted_sampling3():
    """Test sampling with normed abundances but no tail probabilities
    from data."""
    # Abundances are normalized across rule instances involving the same gene.
    os.environ['TEST_FLAG'] = 'TRUE'
    map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    mapk1 = Agent('MAPK1', db_refs={'HGNC': '6871'})
    mapk3 = Agent('MAPK3', db_refs={'HGNC': '6877'})
    jun = Agent('JUN', db_refs={'HGNC': '6204'})
    mapk1_p218 = Agent('MAPK1',
                       mods=[ModCondition('phosphorylation', 'S', '218')],
                       db_refs={'HGNC': '6871'})
    mapk1_p222 = Agent('MAPK1',
                       mods=[ModCondition('phosphorylation', 'S', '222')],
                       db_refs={'HGNC': '6871'})
    mapk3_phos = Agent('MAPK3',
                       mods=[ModCondition('phosphorylation')],
                       db_refs={'HGNC': '6877'})
    st1 = Phosphorylation(map2k1, mapk3)
    st2 = Phosphorylation(map2k1, mapk1, 'S', '218')
    st3 = Phosphorylation(map2k1, mapk1, 'S', '222')
    st4 = Phosphorylation(mapk3_phos, jun)
    st5 = Phosphorylation(mapk1_p218, jun)
    st6 = Phosphorylation(mapk1_p222, jun)
    stmt_to_check = Phosphorylation(map2k1, jun)
    # Make model
    pa = PysbAssembler()
    pa.add_statements([st1, st2, st3, st4, st5, st6])
    pa.make_model(policies='one_step')
    # Set the initial conditions
    mapk1_monomer = pa.model.all_components()['MAPK1']
    mapk3_monomer = pa.model.all_components()['MAPK3']
    set_base_initial_condition(pa.model, mapk1_monomer, 50)
    set_base_initial_condition(pa.model, mapk3_monomer, 50)
    # Do sampling
    mc = PysbModelChecker(pa.model, [stmt_to_check], do_sampling=True, seed=1)
    mc.prune_influence_map()
    results = mc.check_model(max_path_length=5, max_paths=100)
    assert type(results) == list
    assert len(results) == 1
    stmt_tuple = results[0]
    assert len(stmt_tuple) == 2
    assert stmt_tuple[0] == stmt_to_check
    path_result = stmt_tuple[1]
    assert type(path_result) == PathResult
    path_lengths = [len(p) for p in path_result.paths]
    assert max(path_lengths) <= 5
    # There are two distinct paths
    path_ctr = Counter(path_result.paths)
    assert len(path_ctr) == 3
    assert path_ctr[(('MAP2K1_phosphorylation_MAPK3_phospho', 1),
                     ('MAPK3_phospho_phosphorylation_JUN_phospho', 1),
                     ('JUN_phospho_p_obs', 1))] == 49, path_ctr
    assert path_ctr[(('MAP2K1_phosphorylation_MAPK1_S218', 1),
                     ('MAPK1_phosphoS218_phosphorylation_JUN_phospho', 1),
                     ('JUN_phospho_p_obs', 1))] == 31, path_ctr
    assert path_ctr[(('MAP2K1_phosphorylation_MAPK1_S222', 1),
                     ('MAPK1_phosphoS222_phosphorylation_JUN_phospho', 1),
                     ('JUN_phospho_p_obs', 1))] == 20, path_ctr


def test_amount_vs_activation():
    p53 = Agent('TP53', db_refs={'HGNC': '11998'})
    pten1 = Agent('PTEN', bound_conditions=[BoundCondition(p53, True)],
                  db_refs={'HGNC': '9588'})
    pten2 = Agent('PTEN', mods=[ModCondition('modification')],
                  db_refs={'HGNC': '9588'})
    mdm2 = Agent('MDM2', db_refs={'HGNC': '6973'})
    test_stmt = Activation(pten1, p53)
    model_stmts = [IncreaseAmount(pten2, p53),
                   Inhibition(mdm2, p53)]
    # Make model
    pa = PysbAssembler()
    pa.add_statements(model_stmts)
    pa.make_model(policies='one_step')
    # Do sampling
    mc = PysbModelChecker(pa.model, [test_stmt])
    mc.prune_influence_map()
    mc.draw_im('test.pdf')
    results = mc.check_model(max_path_length=1, max_paths=1)
    assert results[0][1].result_code == 'NO_PATHS_FOUND', results


# Test other ModelChecker types
st1 = Activation(Agent('A', db_refs={'HGNC': '1'}),
                 Agent('B', db_refs={'HGNC': '2'}))
st2 = Inhibition(Agent('B', db_refs={'HGNC': '2'}),
                 Agent('D', db_refs={'HGNC': '4'}))
st3 = IncreaseAmount(Agent('C', db_refs={'HGNC': '3'}),
                     Agent('B', db_refs={'HGNC': '2'}))
st4 = DecreaseAmount(Agent('C', db_refs={'HGNC': '3'}),
                     Agent('D', db_refs={'HGNC': '4'}))
st5 = IncreaseAmount(Agent('D', db_refs={'HGNC': '4'}),
                     Agent('E', db_refs={'HGNC': '5'}))
st6 = Inhibition(Agent('A', db_refs={'HGNC': '1'}),
                 Agent('B', db_refs={'HGNC': '2'}))
st7 = DecreaseAmount(Agent('B', db_refs={'HGNC': '2'}),
                     Agent('D', db_refs={'HGNC': '4'}))
st8 = IncreaseAmount(Agent('E', db_refs={'HGNC': '5'}),
                     Agent('B', db_refs={'HGNC': '2'}))
statements = [st1, st2, st3, st4, st5, st6, st7, st8]

test_st1 = Activation(Agent('A', db_refs={'HGNC': '1'}),
                      Agent('E', db_refs={'HGNC': '5'}))
test_st2 = Inhibition(Agent('A', db_refs={'HGNC': '1'}),
                      Agent('E', db_refs={'HGNC': '5'}))
test_st3 = Activation(Agent('A', db_refs={'HGNC': '1'}),
                      Agent('C', db_refs={'HGNC': '3'}))
test_st4 = Activation(Agent('F', db_refs={'HGNC': '6'}),
                      Agent('B', db_refs={'HGNC': '2'}))
test_st5 = DecreaseAmount(Agent('B', db_refs={'HGNC': '2'}),
                          Agent('F', db_refs={'HGNC': '6'}))
test_st6 = ActiveForm(Agent('A', db_refs={'HGNC': '1'}), None, True)
test_st7 = DecreaseAmount(Agent('B', db_refs={'HGNC': '2'}),
                          Agent('B', db_refs={'HGNC': '2'}))
test_statements = [
    test_st1, test_st2, test_st3, test_st4, test_st5, test_st6, test_st7]


def test_unsigned_path():
    ia = IndraNetAssembler(statements)
    unsigned_model = ia.make_model(graph_type='digraph')
    umc = UnsignedGraphModelChecker(unsigned_model)
    umc.add_statements(test_statements)
    results = umc.check_model()
    # Paths found
    assert results[0][1].result_code == 'PATHS_FOUND'
    assert results[0][1].paths[0] == (('A', 0), ('B', 0), ('D', 0), ('E', 0))
    assert results[1][1].result_code == 'PATHS_FOUND'
    assert results[1][1].paths[0] == (('A', 0), ('B', 0), ('D', 0), ('E', 0))
    # Fail cases
    assert results[2][1].result_code == 'NO_PATHS_FOUND'
    assert results[3][1].result_code == 'SUBJECT_NOT_FOUND'
    assert results[4][1].result_code == 'OBJECT_NOT_FOUND'
    assert results[5][1].result_code == 'STATEMENT_TYPE_NOT_HANDLED'
    # Loop paths
    assert results[6][1].result_code == 'PATHS_FOUND'
    assert results[6][1].paths[0] == (('B', 0), ('D', 0), ('E', 0), ('B', 0))
    # Test reporting
    path0 = results[0][1].paths[0]
    path1 = results[1][1].paths[0]
    path6 = results[6][1].paths[0]
    stmts0 = stmts_from_indranet_path(
        path0, unsigned_model, False, False, statements)
    stmts1 = stmts_from_indranet_path(
        path0, unsigned_model, False, False, statements)
    stmts6 = stmts_from_indranet_path(
        path6, unsigned_model, False, False, statements)
    assert stmts0 == stmts1 == [[st1, st6], [st2, st7], [st5]]
    assert stmts6 == [[st2, st7], [st5], [st8]]


def test_signed_path():
    ia = IndraNetAssembler(statements)
    signed_model = ia.make_model(graph_type='signed')
    smc = SignedGraphModelChecker(signed_model)
    smc.add_statements(test_statements)
    results = smc.check_model()
    # Paths found
    assert results[0][1].result_code == 'PATHS_FOUND'
    assert results[0][1].paths[0] == (('A', 0), ('B', 1), ('D', 0), ('E', 0))
    assert results[1][1].result_code == 'PATHS_FOUND'
    assert results[1][1].paths[0] == (('A', 0), ('B', 0), ('D', 1), ('E', 1))
    # Fail cases
    assert results[2][1].result_code == 'NO_PATHS_FOUND'
    assert results[3][1].result_code == 'SUBJECT_NOT_FOUND'
    assert results[4][1].result_code == 'OBJECT_NOT_FOUND'
    assert results[5][1].result_code == 'STATEMENT_TYPE_NOT_HANDLED'
    # Loop paths
    assert results[6][1].result_code == 'PATHS_FOUND'
    assert results[6][1].paths[0] == (('B', 0), ('D', 1), ('E', 1), ('B', 1))
    # Test reporting
    path0 = results[0][1].paths[0]
    path1 = results[1][1].paths[0]
    path6 = results[6][1].paths[0]
    stmts0 = stmts_from_indranet_path(
        path0, signed_model, True, False, statements)
    assert stmts0 == [[st6], [st2, st7], [st5]]
    stmts1 = stmts_from_indranet_path(
        path1, signed_model, True, False, statements)
    assert stmts1 == [[st1], [st2, st7], [st5]]
    stmts6 = stmts_from_indranet_path(
        path6, signed_model, True, False, statements)
    assert stmts6 == [[st2, st7], [st5], [st8]]


def test_pybel_path():
    pba = PybelAssembler(statements)
    pybel_model = pba.make_model()
    pbmc = PybelModelChecker(pybel_model)
    pbmc.add_statements(test_statements)
    results = pbmc.check_model()
    a = _get_agent_node(Agent('A', db_refs={'HGNC': '1'}))[0]
    b = _get_agent_node(Agent('B', db_refs={'HGNC': '2'}))[0]
    d = _get_agent_node(Agent('D', db_refs={'HGNC': '4'}))[0]
    e = _get_agent_node(Agent('E', db_refs={'HGNC': '5'}))[0]
    # Paths found
    assert results[0][1].result_code == 'PATHS_FOUND'
    assert results[0][1].paths[0] == ((a, 0), (b, 1), (d, 0), (e, 0))
    assert results[1][1].result_code == 'PATHS_FOUND'
    assert results[1][1].paths[0] == ((a, 0), (b, 0), (d, 1), (e, 1))
    # Fail cases
    assert results[2][1].result_code == 'NO_PATHS_FOUND'
    assert results[3][1].result_code == 'SUBJECT_NOT_FOUND', results[3]
    assert results[4][1].result_code == 'OBJECT_NOT_FOUND'
    assert results[5][1].result_code == 'STATEMENT_TYPE_NOT_HANDLED'
    # Loop paths
    assert results[6][1].result_code == 'PATHS_FOUND'
    assert results[6][1].paths[0] == ((b, 0), (d, 1), (e, 1), (b, 1))
    # Test reporting
    path0 = results[0][1].paths[0]
    path1 = results[1][1].paths[0]
    path6 = results[6][1].paths[0]
    stmts0 = stmts_from_pybel_path(path0, pybel_model, False, statements)
    assert stmts0 == [[st6], [st2, st7], [st5]], stmts0
    stmts1 = stmts_from_pybel_path(path1, pybel_model, False, statements)
    assert stmts1 == [[st1], [st2, st7], [st5]], stmts1
    stmts6 = stmts_from_pybel_path(path6, pybel_model, False, statements)
    assert stmts6 == [[st2, st7], [st5], [st8]]


def test_pybel_active_form_path():
    braf = Agent('BRAF', db_refs={'HGNC': '1097'})
    map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    map2k1_phos = Agent('MAP2K1', db_refs={'HGNC': '6840'}, mods=[
        ModCondition('phosphorylation', 'S', '218', True),
        ModCondition('phosphorylation', 'S', '222', True)])
    mapk1 = Agent('MAPK1', db_refs={'UP': 'P28482'})
    mapk1_phos = Agent('MAPK1', db_refs={'UP': 'P28482'}, mods=[
        ModCondition('phosphorylation', 'T', '185', True),
        ModCondition('phosphorylation', 'Y', '187', True)])
    elk1 = Agent('ELK1', db_refs={'HGNC': '3321'})
    elk1_phos = Agent('ELK1', db_refs={'HGNC': '3321'}, mods=[
        ModCondition('phosphorylation', 'S', '383', True),
        ModCondition('phosphorylation', 'S', '389', True)])
    stmts = [
        Phosphorylation(braf, map2k1, 'S', '222'),
        ActiveForm(map2k1_phos, 'activity', True),
        Phosphorylation(map2k1, mapk1, 'Y', '187'),
        ActiveForm(mapk1_phos, 'activity', True),
        Phosphorylation(mapk1, elk1, 'S', '383'),
        ActiveForm(elk1_phos, 'activity', True)]
    pba = PybelAssembler(stmts)
    pybel_model = pba.make_model()
    pbmc = PybelModelChecker(pybel_model)
    pbmc.add_statements([Activation(braf, elk1)])
    results = pbmc.check_model(1, 10)
    assert results[0][1].path_found, results
    stmts_from_path = stmts_from_pybel_path(
        results[0][1].paths[0], pybel_model, False, stmts)
    assert stmts_from_path == [[stmt] for stmt in stmts]


def test_refinements():
    prkcb = Agent('PRKCB', db_refs={'TEXT': 'PRKCB', 'HGNC': '9395'})
    gsk3b = Agent('GSK3B', db_refs={'TEXT': 'GSK3B', 'HGNC': '4617'})
    map2k1 = Agent('MAP2K1', db_refs={'TEXT': 'MAP2K1', 'HGNC': '6840'})
    mapk1 = Agent('MAPK1', db_refs={'TEXT': 'MAPK1', 'HGNC': '6871'})
    mek = Agent('MEK', db_refs={'TEXT': 'MEK', 'FPLX': 'MEK'})
    erk = Agent('ERK', db_refs={'TEXT': 'ERK', 'FPLX': 'ERK'})
    # Model statements are refined versions of test statements
    model_stmts = [Phosphorylation(prkcb, gsk3b, 'S', '9'),
                   Phosphorylation(map2k1, mapk1)]
    test_stmts = [Phosphorylation(prkcb, gsk3b, 'S'),
                  Phosphorylation(mek, erk)]
    # Test PyBEL
    pba = PybelAssembler(model_stmts)
    pybel_model = pba.make_model()
    pbmc = PybelModelChecker(pybel_model, test_stmts)
    results = pbmc.check_model()
    assert results[0][1].path_found
    assert results[1][1].path_found
    path0 = results[0][1].paths[0]
    path1 = results[1][1].paths[0]
    assert stmts_from_pybel_path(
        path0, pybel_model, False, model_stmts) == [[model_stmts[0]]]
    assert stmts_from_pybel_path(
        path1, pybel_model, False, model_stmts) == [[model_stmts[1]]]
    # Test PySB
    pa = PysbAssembler(model_stmts)
    pysb_model = pa.make_model()
    pmc = PysbModelChecker(pysb_model, test_stmts, model_stmts=model_stmts)
    results = pmc.check_model()
    assert results[0][1].path_found
    assert results[1][1].path_found
    path0 = results[0][1].paths[0]
    path1 = results[1][1].paths[0]
    assert stmts_from_pysb_path(
        path0, pysb_model, model_stmts) == [model_stmts[0]]
    assert stmts_from_pysb_path(
        path1, pysb_model, model_stmts) == [model_stmts[1]]
    # Test unsigned graph
    ia = IndraNetAssembler(model_stmts)
    unsigned_model = ia.make_model(graph_type='digraph')
    umc = UnsignedGraphModelChecker(unsigned_model, test_stmts)
    results = umc.check_model()
    assert results[0][1].path_found
    assert results[1][1].path_found
    path0 = results[0][1].paths[0]
    path1 = results[1][1].paths[0]
    assert stmts_from_indranet_path(
        path0, unsigned_model, False, False, model_stmts) == [[model_stmts[0]]]
    assert stmts_from_indranet_path(
        path1, unsigned_model, False, False, model_stmts) == [[model_stmts[1]]]
    # Test signed graph
    # Can't use phosphorylations for signed graph, using activations
    model_stmts = [Activation(map2k1, mapk1)]
    test_stmts = [Activation(mek, erk)]
    ia = IndraNetAssembler(model_stmts)
    signed_model = ia.make_model(graph_type='signed')
    smc = SignedGraphModelChecker(signed_model, test_stmts)
    results = smc.check_model()
    assert results[0][1].path_found
    path0 = results[0][1].paths[0]
    assert stmts_from_indranet_path(
        path0, signed_model, True, False, model_stmts) == [[model_stmts[0]]]


def test_pybel_edge_types():
    a = Agent('A', db_refs={'HGNC': '1'})
    b = Agent('B', db_refs={'HGNC': '2'})
    c = Agent('C', db_refs={'HGNC': '3'})
    a_b = Agent('A', db_refs={'HGNC': '1'}, bound_conditions=[BoundCondition(b)])
    c_phos = Agent('C', db_refs={'HGNC': '3'},
                   mods=[ModCondition('phosphorylation', 'S', '218')])
    model_stmts = [Complex([a, b]),
                   Phosphorylation(b, c, 'S', '218')]
    test_stmts = [Activation(a, c),
                  Activation(b, c),
                  Activation(b, a),
                  Activation(b, a_b),
                  Activation(c, c_phos)]
    pba = PybelAssembler(model_stmts)
    pybel_model = pba.make_model()
    pbmc = PybelModelChecker(pybel_model, test_stmts)
    # Do not include hasVariant and partOf edges at all
    pbmc.graph = None
    pbmc.get_graph(include_variants=False, symmetric_variant_links=False,
                   include_components=False, symmetric_component_links=False)
    results = pbmc.check_model()
    assert not results[0][1].path_found
    assert not results[1][1].path_found
    assert not results[2][1].path_found
    assert not results[3][1].path_found, results[3][1]
    assert not results[4][1].path_found
    # Include hasVariant and partOf edges without symmetric links
    pbmc.graph = None
    pbmc.get_graph(include_variants=True, symmetric_variant_links=False,
                   include_components=True, symmetric_component_links=False)
    results = pbmc.check_model()
    assert not results[0][1].path_found
    assert not results[1][1].path_found
    assert not results[2][1].path_found
    assert results[3][1].path_found
    assert results[4][1].path_found
    path3 = results[3][1].paths[0]
    path4 = results[4][1].paths[0]
    path_stmt_3 = stmts_from_pybel_path(
        path3, pybel_model, False, model_stmts)[0][0]
    path_stmt_4 = stmts_from_pybel_path(
        path4, pybel_model, False, model_stmts)[0][0]
    assert isinstance(path_stmt_3, PybelEdge)
    assert path_stmt_3.relation == 'partOf'
    assert isinstance(path_stmt_4, PybelEdge)
    assert path_stmt_4.relation == 'hasVariant'
    # Include symmetric links
    pbmc.graph = None
    pbmc.get_graph(include_variants=True, symmetric_variant_links=True,
                   include_components=True, symmetric_component_links=True)
    results = pbmc.check_model()
    assert results[0][1].path_found
    assert results[1][1].path_found
    assert results[2][1].path_found
    assert results[3][1].path_found
    assert results[4][1].path_found


def test_pybel_edge_to_english():
    pe = PybelEdge(
        Agent('EGF', bound_conditions=[BoundCondition(Agent('EGFR'))]),
        Agent('EGF'), 'partOf', True)
    s = pybel_edge_to_english(pe)
    assert s == 'EGF bound to EGFR has a component EGF.'
    pe = PybelEdge(
        Agent('EGF'),
        Agent('EGF', bound_conditions=[BoundCondition(Agent('EGFR'))]),
        'partOf', False)
    s = pybel_edge_to_english(pe)
    assert s == 'EGF is a part of EGF bound to EGFR.'
    pe = PybelEdge(
        Agent('BRAF'),
        Agent('BRAF', mods=[ModCondition('phosphorylation', 'T', '396')]),
        'hasVariant', False)
    s = pybel_edge_to_english(pe)
    assert s == 'BRAF has a variant BRAF phosphorylated on T396.'
    pe = PybelEdge(
        Agent('BRAF', mods=[ModCondition('phosphorylation', 'T', '396')]),
        Agent('BRAF'),
        'hasVariant', True)
    s = pybel_edge_to_english(pe)
    assert s == 'BRAF phosphorylated on T396 is a variant of BRAF.'


# Test graph conversion
def test_signed_edges_to_nodes():
    edge_dict = {'extra_data': {'list': ['value'], 'float': 0.123456},
                 'weight': 0.987654}
    g = nx.MultiDiGraph()
    g.add_edge('A', 'B', sign=0, **edge_dict)
    g.add_edge('A', 'C', sign=0, **edge_dict)
    g.add_edge('B', 'D', sign=1, **edge_dict)
    g.add_edge('C', 'D', sign=0, **edge_dict)
    g.add_edge('E', 'B', sign=1, **edge_dict)
    assert len(g.edges) == 5
    assert len(g.nodes) == 5
    # Create a signed nodes graph without pruning
    sng = signed_edges_to_signed_nodes(g, prune_nodes=False)
    assert len(sng.nodes) == 10
    assert len(sng.edges) == 10
    assert set(sng.nodes) == {('A', 0), ('A', 1), ('B', 0), ('B', 1), ('C', 0),
                              ('C', 1), ('D', 0), ('D', 1), ('E', 0), ('E', 1)}
    assert (('A', 0), ('B', 0)) in sng.edges
    assert (('A', 1), ('B', 1)) in sng.edges
    assert (('B', 0), ('D', 1)) in sng.edges
    assert (('B', 1), ('D', 0)) in sng.edges
    # Create a signed nodes graph with pruning
    psng = signed_edges_to_signed_nodes(g, prune_nodes=True)
    assert len(psng.nodes) == 7
    assert ('A', 1) not in psng.nodes
    assert ('C', 1) not in psng.nodes
    assert ('E', 1) not in psng.nodes
    assert len(psng.edges) == 6
    # Create a signed nodes graph with weight data with pruning
    psng_wed = signed_edges_to_signed_nodes(g, prune_nodes=True,
                                            copy_edge_data={'weight'})
    for edge in psng_wed.edges:
        assert psng_wed.edges[edge]['weight'] == 0.987654
        assert 'sign' not in psng_wed.edges[edge]
        assert 'extra_data' not in psng_wed.edges[edge]

    # Create a signed nodes graph with all edge data with pruning
    psng_ed = signed_edges_to_signed_nodes(g, prune_nodes=True,
                                            copy_edge_data=True)
    for edge in psng_ed.edges:
        assert psng_ed.edges[edge]['weight'] == 0.987654
        assert 'sign' not in psng_ed.edges[edge]
        assert 'extra_data' in psng_ed.edges[edge],\
            psng_ed.edges[edge].items()
        assert psng_ed.edges[edge]['extra_data']['list'] == ['value']
        assert psng_ed.edges[edge]['extra_data']['float'] == 0.123456


def test_path_fixed_length():
    model_stmts = [
        IncreaseAmount(Agent('A', db_refs={'HGNC': '1'}),
                       Agent('B', db_refs={'HGNC': '2'})),
        IncreaseAmount(Agent('B', db_refs={'HGNC': '2'}),
                       Agent('C', db_refs={'HGNC': '3'})),
        IncreaseAmount(Agent('C', db_refs={'HGNC': '3'}),
                       Agent('D', db_refs={'HGNC': '4'})),
        IncreaseAmount(Agent('A', db_refs={'HGNC': '1'}),
                       Agent('C', db_refs={'HGNC': '3'})),
        IncreaseAmount(Agent('B', db_refs={'HGNC': '2'}),
                       Agent('D', db_refs={'HGNC': '4'}))
    ]
    test_stmt = IncreaseAmount(Agent('A', db_refs={'HGNC': 1}),
                               Agent('D', db_refs={'HGNC': 4}))
    # There are two two-step paths and one three-step path.
    # ModelChecker should return all two-step paths.
    ia = IndraNetAssembler(model_stmts)
    unsigned_model = ia.make_model(graph_type='digraph')
    umc = UnsignedGraphModelChecker(unsigned_model)
    res = umc.check_statement(test_stmt, max_paths=10, max_path_length=5)
    assert res.path_found
    assert len(res.paths) == 2, len(res.paths)
    assert len(res.paths[0]) == len(res.paths[1]) == 3  # 3 nodes = 2 edges/steps


if __name__ == '__main__':
    test_prune_influence_map_subj_obj()

# TODO Add tests for autophosphorylation
# TODO Add test for transphosphorylation

# FIXME: Issue: Increasing kinase activity doesn't make it capable of executing
# phosphorylation statements
# FIXME Issue increase activity (generic) doesn't make something capable of
# executing phospho (or other statements)


# Goal: Be able to check generic phosphorylations against specific rules
# and vice versa.
# 0. (('T185', 'p'),('Y187','p'), ...,), 'kinase', 'is_active')
# 1. (('T185', 'p'),('Y187','p'), ...,), 'kinase', 'is_inactive')
# 2. How to figure out that 'kinase'
# 1. Add activity states to agents, and have check_activation work with
#    grounded monomers
# 2. Also need activity observable matching to account for different types of
#    activity, so that different types of activity can be checked
# 2. Need to be able to annotate specific site/state combinations as
#    active forms
# 3. Need to make grounded_mp generation work with MutConditions and
#    bound conditions (bound conditions would need to check for bonds to
# 4. Grounded monomer patterns for check_activation
# 5. Check for complexes (using contact map? Or express the complex as an
#    observable and look for paths?
# Issues--if a rule activity isn't contingent on a particular mod,
# then were will be no causal connection between any upstream modifying
# rules and the downstream rule.

# Active Form vs. Conditions on enzymes
# Identify conflicting ModConditions?
# Refactor out active forms as a preassembly step

# Add "active" into assembled Pysb models as a requirement of every enzyme,
# and add ActiveForms in as rules in the model. This has the advantage of
# working when the activeform is not known (allows tying together activation
# and mods). Problem is what to do with the rules that have mod conditions
# already--add active in? Active flag approach has the advantage of a single
# role for a each substrate, which (I think) would prevent an explosion of
# paths

# Simplest thing: take any rule where the enzyme has no conditions on it
# and replace it with a comparable rule for each activeform.

# Create rules for each active form

# So what's needed is an assembly procedure where an active form is applied
# across all rules for which that protein is the enzyme.

# Need to know which agent is the "enzyme" in rules, so that these can
# be prioritized, and activeforms applied

# Apply weights based on evidence/belief scores;
# Apply weights based on positive

# Another issue--doesn't know that RAF1(phospho='p') should be satisfied
# by RAF1(S259='p'). A big problem, even after pre-assembly--because longer
# paths where a 'phospho' is the observable will never be satisfied.

# Does this mean that we need a PySB ComplexPattern -> Agent mapping, that
# we can subsequently use for refinements?

# Need to handle complex statements. Would show that one_step approach
# would not satisfy constraint, but two-step approach could, where the
# Complex information was specified.
# Can probably handle all modifications in a generic function.
# Then need to handle: Complex, Dephosphorylation.
# Then Gef/Gap?
# Then Activation/ActiveForm.
# Get the stuff from databases involving the canonical proteins,
# and show that a simple model satisfies it.
# Try to build the model using natural language?
#
# By tying molecules to biological processes, we can even check that
# these types of high-level observations are satisfied.
#
# Need to handle case where Phosphorylation site is not specified by
# statement, but is actually handled in the model (i.e., need to know
# that a particular site name and state corresponds to a phosphorylation.
# Points to need to have an additional data structure annotating agents,
# sites, states.
#
# Need to handle embeddings of complex patterns where sites can have both
# modification state and bonds
#
# Need to handle reversible rules!
#
# Should probably build in some way of returning the paths found
#
# Save all the paths that a particular rule is on--then if you're wondering
# why it's in the model, you look at all of the statements for which that
# rule provides a path  .
#
# When Ras machine finds a new finding, it can be checked to see if it's
# satisfied by the model.
