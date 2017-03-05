from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers import PysbAssembler
from indra.assemblers import pysb_assembler as pa
from indra.assemblers.pysb_assembler import PysbPreassembler
from indra.statements import *
from pysb import bng, WILD, Monomer, Annotation
from pysb.testing import with_model

def test_pysb_assembler_complex1():
    member1 = Agent('BRAF')
    member2 = Agent('MEK1')
    stmt = Complex([member1, member2])
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==2)
    assert(len(model.monomers)==2)

def test_pysb_assembler_complex2():
    member1 = Agent('BRAF')
    member2 = Agent('MEK1')
    member3 = Agent('ERK1')
    stmt = Complex([member1, member2, member3])
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==6)
    assert(len(model.monomers)==3)

def test_pysb_assembler_complex3():
    hras = Agent('HRAS')
    member1 = Agent('BRAF', bound_conditions=[BoundCondition(hras, True)])
    member2 = Agent('MEK1')
    stmt = Complex([member1, member2])
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==2)
    assert(len(model.monomers)==3)

def test_pysb_assembler_complex_twostep():
    member1 = Agent('BRAF')
    member2 = Agent('MEK1')
    stmt = Complex([member1, member2])
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model(policies='two_step')
    assert(len(model.rules)==2)
    assert(len(model.monomers)==2)

def test_pysb_assembler_complex_multiway():
    member1 = Agent('BRAF')
    member2 = Agent('MEK1')
    member3 = Agent('ERK1')
    stmt = Complex([member1, member2, member3])
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model(policies='multi_way')
    assert(len(model.rules)==2)
    assert(len(model.monomers)==3)

def test_pysb_assembler_actsub():
    stmt = ActiveForm(Agent('BRAF', mutations=[MutCondition('600', 'V', 'E')]),
                      'activity', True)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model(policies='two_step')
    assert(len(model.rules)==0)
    assert(len(model.monomers)==1)

def test_pysb_assembler_phos_noenz():
    enz = None
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==0)
    assert(len(model.monomers)==0)

def test_pysb_assembler_dephos_noenz():
    enz = None
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==0)
    assert(len(model.monomers)==0)

def test_pysb_assembler_phos1():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_phos2():
    hras = Agent('HRAS')
    enz = Agent('BRAF', bound_conditions=[BoundCondition(hras, True)])
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==3)

def test_pysb_assembler_phos3():
    hras = Agent('HRAS')
    erk1 = Agent('ERK1')
    enz = Agent('BRAF', bound_conditions=[BoundCondition(hras, True)])
    sub = Agent('MEK1', bound_conditions=[BoundCondition(erk1, True)])
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==4)

def test_pysb_assembler_phos4():
    hras = Agent('HRAS')
    erk1 = Agent('ERK1')
    enz = Agent('BRAF', bound_conditions=[BoundCondition(hras, True)])
    sub = Agent('MEK1', bound_conditions=[BoundCondition(erk1, False)])
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==4)

def test_pysb_assembler_autophos1():
    enz = Agent('MEK1')
    stmt = Autophosphorylation(enz, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_autophos2():
    raf1 = Agent('RAF1')
    enz = Agent('MEK1', bound_conditions=[BoundCondition(raf1, True)])
    stmt = Autophosphorylation(enz, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_autophos3():
    egfr = Agent('EGFR')
    enz = Agent('EGFR', bound_conditions=[BoundCondition(egfr, True)])
    stmt = Autophosphorylation(enz, 'tyrosine')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_transphos1():
    egfr = Agent('EGFR')
    enz = Agent('EGFR', bound_conditions=[BoundCondition(egfr, True)])
    stmt = Transphosphorylation(enz, 'tyrosine')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_act1():
    egfr = Agent('EGFR')
    subj = Agent('GRB2', bound_conditions=[BoundCondition(egfr, True)])
    obj = Agent('SOS1')
    stmt = Activation(subj, obj)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==3)

def test_pysb_assembler_dephos1():
    phos = Agent('PP2A')
    sub = Agent('MEK1')
    stmt = Dephosphorylation(phos, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_dephos2():
    phos = Agent('PP2A')
    raf1 = Agent('RAF1')
    sub = Agent('MEK1', bound_conditions=[BoundCondition(raf1, True)])
    stmt = Dephosphorylation(phos, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==3)

def test_pysb_assembler_rasgef1():
    gef = Agent('SOS1')
    ras = Agent('HRAS')
    stmt = RasGef(gef, ras)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_rasgap1():
    gap = Agent('NF1')
    ras = Agent('HRAS')
    stmt = RasGap(gap, ras)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_actmod1():
    mek = Agent('MEK')
    erk = Agent('ERK')
    stmts = []
    mc1 = ModCondition('phosphorylation', 'serine', '218')
    mc2 = ModCondition('phosphorylation', 'serine', '222')
    stmts.append(ActiveForm(Agent('MEK', mods=[mc1, mc2]), 'activity', True))
    stmts.append(Phosphorylation(mek, erk, 'threonine', '185'))
    stmts.append(Phosphorylation(mek, erk, 'tyrosine', '187'))
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    assert(len(model.rules)==2)
    assert(len(model.monomers)==2)
    model = pa.make_model(policies='two_step')
    assert(len(model.rules)==5)

def test_pysb_assembler_actmod2():
    mek = Agent('MEK', activity=ActivityCondition('activity', True))
    erk = Agent('ERK')
    stmts = []
    stmts.append(ActiveForm(Agent('MEK',
                    mods=[ModCondition('phosphorylation', 'serine', '218')]),
                    'activity', True))
    stmts.append(ActiveForm(Agent('MEK',
                    mods=[ModCondition('phosphorylation', 'serine', '222')]),
                    'activity', True))
    stmts.append(Phosphorylation(mek, erk, 'threonine', '185'))
    stmts.append(Phosphorylation(mek, erk, 'tyrosine', '187'))
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    assert(len(model.rules)==4)
    assert(len(model.monomers)==2)
    model = pa.make_model(policies='two_step')
    assert(len(model.rules)==9)

def test_pysb_assembler_phos_twostep1():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler(policies='two_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==3)
    assert(len(model.monomers)==2)

def test_pysb_assembler_twostep_mixed():
    member1 = Agent('BRAF')
    member2 = Agent('RAF1')
    st1 = Complex([member1, member2])
    st2 = Phosphorylation(Agent('MAP2K1'), Agent('MAPK3'))
    pa = PysbAssembler()
    pa.add_statements([st1, st2])
    pa.make_model(policies='two_step')
    assert(len(pa.model.rules)==5)
    assert(len(pa.model.monomers)==4)

def test_pysb_assembler_phos_twostep_local():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model(policies='two_step')
    assert(len(model.rules)==3)
    assert(len(model.monomers)==2)

def test_pysb_assembler_phos_twostep_local_to_global():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'serine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model(policies='two_step')
    # This call should have reverted to default policy
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_dephos_twostep1():
    phos = Agent('PP2A')
    sub = Agent('MEK1')
    stmt = Dephosphorylation(phos, sub, 'serine', '222')
    pa = PysbAssembler(policies='two_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==3)
    assert(len(model.monomers)==2)

def test_statement_specific_policies():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    phos = Agent('PP2A')
    stmt1 = Phosphorylation(enz, sub, 'serine', '222')
    stmt2 = Dephosphorylation(phos, sub, 'serine', '222')
    policies = {'Phosphorylation': 'two_step',
                'Dephosphorylation': 'interactions_only'}
    pa = PysbAssembler(policies=policies)
    pa.add_statements([stmt1, stmt2])
    model = pa.make_model()
    assert(len(model.rules)==4)
    assert(len(model.monomers)==3)

def test_unspecified_statement_policies():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    phos = Agent('PP2A')
    stmt1 = Phosphorylation(enz, sub, 'serine', '222')
    stmt2 = Dephosphorylation(phos, sub, 'serine', '222')
    policies = {'Phosphorylation': 'two_step',
                'other': 'interactions_only'}
    pa = PysbAssembler(policies=policies)
    pa.add_statements([stmt1, stmt2])
    model = pa.make_model()
    assert(len(model.rules)==4)
    assert(len(model.monomers)==3)

def test_activity_activity():
    subj = Agent('KRAS')
    obj = Agent('BRAF')
    stmt = Activation(subj, obj)
    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_activity_activity2():
    subj = Agent('KRAS')
    obj = Agent('BRAF')
    stmt = Activation(subj, obj)
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_activity_activity2():
    subj = Agent('Vemurafenib')
    obj = Agent('BRAF')
    stmt = Inhibition(subj, obj)
    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_activity_activity3():
    subj = Agent('Vemurafenib')
    obj = Agent('BRAF')
    stmt = Inhibition(subj, obj)
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_rule_name_str_1():
    s = pa.get_agent_rule_str(Agent('BRAF'))
    assert(s == 'BRAF')

def test_rule_name_str_2():
    a = Agent('GRB2',
              bound_conditions=[BoundCondition(Agent('EGFR'), True)])
    s = pa.get_agent_rule_str(a)
    assert(s == 'GRB2_EGFR')

def test_rule_name_str_3():
    a = Agent('GRB2',
              bound_conditions=[BoundCondition(Agent('EGFR'), False)])
    s = pa.get_agent_rule_str(a)
    assert(s == 'GRB2_nEGFR')

def test_rule_name_str_4():
    a = Agent('BRAF', mods=[ModCondition('phosphorylation', 'serine')])
    s = pa.get_agent_rule_str(a)
    assert(s == 'BRAF_phosphoS')

def test_rule_name_str_5():
    a = Agent('BRAF', mods=[ModCondition('phosphorylation', 'serine', '123')])
    s = pa.get_agent_rule_str(a)
    assert(s == 'BRAF_phosphoS123')

def test_neg_act_mod():
    mc = ModCondition('phosphorylation', 'serine', '123', False)
    st1 = ActiveForm(Agent('BRAF', mods=[mc]), 'activity', True)
    braf = Agent('BRAF', activity=ActivityCondition('active', True))
    st2 = Phosphorylation(braf, Agent('MAP2K2'))
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([st1, st2])
    pa.make_model()
    assert(len(pa.model.rules) == 1)
    r = pa.model.rules[0]
    braf = r.reactant_pattern.complex_patterns[0].monomer_patterns[0]
    assert(braf.monomer.name == 'BRAF')
    assert(braf.site_conditions == {'S123': ('u', WILD)})

def test_pos_agent_mod():
    mc = ModCondition('phosphorylation', 'serine', '123', True)
    st = Phosphorylation(Agent('BRAF', mods=[mc]), Agent('MAP2K2'))
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([st])
    pa.make_model()
    assert(len(pa.model.rules) == 1)
    r = pa.model.rules[0]
    braf = r.reactant_pattern.complex_patterns[0].monomer_patterns[0]
    assert(braf.monomer.name == 'BRAF')
    assert(braf.site_conditions == {'S123': ('p', WILD)})

def test_neg_agent_mod():
    mc = ModCondition('phosphorylation', 'serine', '123', False)
    st = Phosphorylation(Agent('BRAF', mods=[mc]), Agent('MAP2K2'))
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([st])
    pa.make_model()
    assert(len(pa.model.rules) == 1)
    r = pa.model.rules[0]
    braf = r.reactant_pattern.complex_patterns[0].monomer_patterns[0]
    assert(braf.monomer.name == 'BRAF')
    assert(braf.site_conditions == {'S123': ('u', WILD)})

def test_mut():
    mut = MutCondition('600', 'V', 'E')
    st = Phosphorylation(Agent('BRAF', mutations=[mut]), Agent('MEK'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(len(pa.model.rules) == 1)
    r = pa.model.rules[0]
    braf = r.reactant_pattern.complex_patterns[0].monomer_patterns[0]
    assert(braf.monomer.name == 'BRAF')
    assert(braf.site_conditions == {'V600': 'E'})

def test_agent_loc():
    st = Phosphorylation(Agent('BRAF', location='cytoplasm'), Agent('MEK'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(len(pa.model.rules) == 1)
    r = pa.model.rules[0]
    braf = r.reactant_pattern.complex_patterns[0].monomer_patterns[0]
    assert(braf.site_conditions == {'loc': 'cytoplasm'})

def test_translocation():
    st = Translocation(Agent('FOXO3A'), 'nucleus', 'cytoplasm')
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(len(pa.model.rules) == 1)
    r = pa.model.rules[0]
    f1 = r.reactant_pattern.complex_patterns[0].monomer_patterns[0]
    assert(f1.site_conditions == {'loc': 'nucleus'})
    f2 = r.product_pattern.complex_patterns[0].monomer_patterns[0]
    assert(f2.site_conditions == {'loc': 'cytoplasm'})
    assert(r.rate_forward.name == 'kf_foxo3a_nucleus_cytoplasm_1')

def test_phos_atpdep():
    st = Phosphorylation(Agent('BRAF'), Agent('MEK'), 'S', '222')
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model(policies='atp_dependent')
    assert(len(pa.model.rules) == 5)

def test_set_context():
    st = Phosphorylation(Agent('MAP2K1'), Agent('MAPK3'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(pa.model.parameters['MAP2K1_0'].value == pa.default_initial_amount)
    assert(pa.model.parameters['MAPK3_0'].value == pa.default_initial_amount)
    pa.set_context('A375_SKIN')
    assert(pa.model.parameters['MAP2K1_0'].value > 10000)
    assert(pa.model.parameters['MAPK3_0'].value > 10000)

def test_set_context_monomer_notfound():
    st = Phosphorylation(Agent('MAP2K1'), Agent('XYZ'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(pa.model.parameters['MAP2K1_0'].value == pa.default_initial_amount)
    assert(pa.model.parameters['XYZ_0'].value == pa.default_initial_amount)
    pa.add_default_initial_conditions(100)
    assert(pa.model.parameters['MAP2K1_0'].value == 100)
    assert(pa.model.parameters['XYZ_0'].value == 100)
    pa.set_context('A375_SKIN')
    assert(pa.model.parameters['MAP2K1_0'].value > 10000)
    assert(pa.model.parameters['XYZ_0'].value == pa.default_initial_amount)

def test_set_context_celltype_notfound():
    st = Phosphorylation(Agent('MAP2K1'), Agent('MAPK3'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    pa.set_context('XYZ')

def test_annotation():
    st = Phosphorylation(Agent('BRAF', db_refs = {'UP': 'P15056'}),
                         Agent('MAP2K2', db_refs = {'HGNC': '6842'}))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(len(pa.model.annotations) == 5)

def test_print_model():
    st = Phosphorylation(Agent('MAP2K1'), Agent('MAPK3'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    pa.save_model('/dev/null')

def test_save_rst():
    st = Phosphorylation(Agent('MAP2K1'), Agent('MAPK3'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    pa.save_rst('/dev/null')

def test_export_model():
    st = Phosphorylation(Agent('MAP2K1'), Agent('MAPK3'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    exp_str = pa.export_model('kappa')
    assert(exp_str)
    exp_str = pa.export_model('bngl')
    assert(exp_str)
    exp_str = pa.export_model('sbml', file_name='/dev/null')
    assert(exp_str)

def test_name_standardize():
    n = pa._n('.*/- ^&#@$')
    assert(isinstance(n, str))
    assert(n == '__________')
    n = pa._n('14-3-3')
    assert(isinstance(n, str))
    assert(n == 'p14_3_3')
    n = pa._n('\U0001F4A9bar')
    assert(isinstance(n, str))
    assert(n == 'bar')

def test_generate_equations():
    st = Phosphorylation(Agent('MAP2K1'), Agent('MAPK3'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    bng.generate_equations(pa.model)

def test_non_python_name_phos():
    st = Phosphorylation(Agent('14-3-3'), Agent('BRAF kinase'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    names = [m.name for m in pa.model.monomers]
    assert('BRAF_kinase' in names)
    assert('p14_3_3' in names)
    bng.generate_equations(pa.model)

def test_non_python_name_bind():
    st = Complex([Agent('14-3-3'), Agent('BRAF kinase')])
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    bng.generate_equations(pa.model)

def test_decreaseamount_one_step():
    subj = Agent('KRAS')
    obj = Agent('BRAF')
    st1 = DecreaseAmount(subj, obj)
    st2 = DecreaseAmount(None, obj)
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([st1, st2])
    model = pa.make_model()
    assert(len(model.rules)==2)
    assert(len(model.monomers)==2)

def test_decreaseamount_interactions_only():
    subj = Agent('KRAS')
    obj = Agent('BRAF')
    st1 = DecreaseAmount(subj, obj)
    st2 = DecreaseAmount(None, obj)
    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([st1, st2])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_increaseamount_one_step():
    subj = Agent('KRAS')
    obj = Agent('BRAF')
    st1 = IncreaseAmount(subj, obj)
    st2 = IncreaseAmount(None, obj)
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([st1, st2])
    model = pa.make_model()
    assert(len(model.rules)==2)
    assert(len(model.monomers)==2)

def test_increaseamount_interactions_only():
    subj = Agent('KRAS')
    obj = Agent('BRAF')
    st1 = IncreaseAmount(subj, obj)
    st2 = IncreaseAmount(None, obj)
    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([st1, st2])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_missing_catalytic_default_site():
    c8 = Agent('CASP8', activity=ActivityCondition('catalytic', True))
    c3 = Agent('CASP3')
    stmt = Activation(c8, c3, 'catalytic')
    # Interactions only
    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([stmt])
    model = pa.make_model()
    # One step
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    # Two step
    pa = PysbAssembler(policies='two_step')
    pa.add_statements([stmt])
    model = pa.make_model()

def test_missing_transcription_default_site():
    p53 = Agent('TP53', activity=ActivityCondition('transcription', True))
    bax = Agent('BAX')
    stmt = Activation(p53, bax)
    # Interactions only
    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([stmt])
    model = pa.make_model()
    # One step
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    # Two step
    pa = PysbAssembler(policies='two_step')
    pa.add_statements([stmt])
    model = pa.make_model()

def test_translocation_loc_special_char():
    st = Translocation(Agent('KSR1'), 'cytoplasm', 'cell surface')
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(len(pa.model.rules) == 1)
    r = pa.model.rules[0]
    f1 = r.reactant_pattern.complex_patterns[0].monomer_patterns[0]
    assert(f1.site_conditions == {'loc': 'cytoplasm'})
    f2 = r.product_pattern.complex_patterns[0].monomer_patterns[0]
    assert(f2.site_conditions == {'loc': 'cell_surface'})
    assert(r.rate_forward.name == 'kf_ksr1_cytoplasm_cell_surface_1')

def test_parse_identifiers_url():
    url1 = 'http://identifiers.org/foo/bar'
    url2 = 'http://identifiers.org/hgnc/12345'
    url3 = 'http://identifiers.org/hgnc/HGNC:12345'
    url4 = 'http://identifiers.org/uniprot/12345'
    url5 = 'http://identifiers.org/chebi/12345'
    url6 = 'http://identifiers.org/interpro/12345'
    url7 = 'http://identifiers.org/pfam/12345'
    (ns, id) = pa.parse_identifiers_url(url1)
    assert ns is None and id is None
    (ns, id) = pa.parse_identifiers_url(url2)
    assert ns is None and id is None
    (ns, id) = pa.parse_identifiers_url(url3)
    assert ns == 'HGNC' and id == '12345'
    (ns, id) = pa.parse_identifiers_url(url4)
    assert ns == 'UP' and id == '12345'
    (ns, id) = pa.parse_identifiers_url(url5)
    assert ns == 'CHEBI' and id == '12345'
    (ns, id) = pa.parse_identifiers_url(url6)
    assert ns == 'IP' and id == '12345'
    (ns, id) = pa.parse_identifiers_url(url7)
    assert ns == 'XFAM' and id == '12345'

@with_model
def test_get_mp_with_grounding():
    foo = Agent('Foo', db_refs={'HGNC': 'foo'})
    a = Agent('A', db_refs={'HGNC': '6840'})
    b = Agent('B', db_refs={'HGNC': '6871'})
    Monomer('A_monomer')
    Monomer('B_monomer')
    Annotation(A_monomer, 'http://identifiers.org/hgnc/HGNC:6840')
    Annotation(B_monomer, 'http://identifiers.org/hgnc/HGNC:6871')
    mps = list(pa.grounded_monomer_patterns(model, foo))
    assert len(mps) == 0
    mps = list(pa.grounded_monomer_patterns(model, a))
    assert len(mps) == 1
    assert mps[0].monomer == A_monomer
    mps = list(pa.grounded_monomer_patterns(model, b))
    assert len(mps) == 1
    assert mps[0].monomer == B_monomer

@with_model
def test_get_mp_with_grounding_2():
    a1 = Agent('A', mods=[ModCondition('phosphorylation', None, None)],
                db_refs={'HGNC': '6840'})
    a2 = Agent('A', mods=[ModCondition('phosphorylation', 'Y', '187')],
                db_refs={'HGNC': '6840'})
    Monomer('A_monomer', ['phospho', 'T185', 'Y187'],
            {'phospho': 'y', 'T185':['u', 'p'], 'Y187':['u','p']})
    Annotation(A_monomer, 'http://identifiers.org/hgnc/HGNC:6840')
    A_monomer.site_annotations = [
        Annotation(('phospho', 'y'), 'phosphorylation', 'is_modification'),
        Annotation(('T185', 'p'), 'phosphorylation', 'is_modification'),
        Annotation(('Y187', 'p'), 'phosphorylation', 'is_modification'),
        Annotation('T185', 'T', 'is_residue'),
        Annotation('T185', '185', 'is_position'),
        Annotation('Y187', 'Y', 'is_residue'),
        Annotation('Y187', '187', 'is_position')
    ]
    mps_1 = list(pa.grounded_monomer_patterns(model, a1))
    assert len(mps_1) == 3
    mps_2 = list(pa.grounded_monomer_patterns(model, a2))
    assert len(mps_2) == 1
    mp = mps_2[0]
    assert mp.monomer == A_monomer
    assert mp.site_conditions == {'Y187': ('p', WILD)}
    # TODO Add test for unmodified agent!
    # TODO Add test involving multiple (possibly degenerate) modifications!
    # TODO Add test for generic double phosphorylation

def test_phospho_assemble_grounding():
    a = Agent('MEK1', db_refs={'HGNC': '6840'})
    b = Agent('ERK2', db_refs={'HGNC': '6871'})
    b_phos = Agent('Foo', mods=[ModCondition('phosphorylation', None, None)],
                    db_refs={'HGNC': '6871'})
    st1 = Phosphorylation(a, b, 'T', '185')
    # One step
    def check_policy(policy):
        pysb_asmb = pa.PysbAssembler(policies=policy)
        pysb_asmb.add_statements([st1])
        model = pysb_asmb.make_model()
        mps = list(pa.grounded_monomer_patterns(model, b_phos))
        assert len(mps) == 1
        assert mps[0].monomer.name == 'ERK2'
        assert mps[0].site_conditions == {'T185': ('p', WILD)}
    for policy in ('one_step', 'interactions_only', 'two_step',
                   'atp_dependent'):
        check_policy(policy)

def test_phospho_mod_grounding():
    a = Agent('MEK1', mods=[ModCondition('phosphorylation', 'S', '218'),
                            ModCondition('phosphorylation', 'S', '222')],
              db_refs={'HGNC': '6840'})
    b = Agent('ERK2', db_refs={'HGNC': '6871'})
    a_phos = Agent('Foo', mods=[ModCondition('phosphorylation', None, None)],
                    db_refs={'HGNC': '6840'})
    st1 = Phosphorylation(a, b, 'T', '185')
    pysb_asmb = pa.PysbAssembler(policies='one_step')
    pysb_asmb.add_statements([st1])
    model = pysb_asmb.make_model()
    mps = list(pa.grounded_monomer_patterns(model, a_phos))
    assert len(mps) == 2
    assert mps[0].monomer.name == 'MEK1'
    assert mps[1].monomer.name == 'MEK1'
    sc = [mp.site_conditions for mp in mps]
    assert {'S218': ('p', WILD)} in sc
    assert {'S222': ('p', WILD)} in sc

def _check_mod_assembly(mod_class):
    subj = Agent('KRAS')
    obj = Agent('BRAF')
    st1 = mod_class(subj, obj)

    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([st1])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

    pa = PysbAssembler(policies='one_step')
    pa.add_statements([st1])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

    pa = PysbAssembler(policies='two_step')
    pa.add_statements([st1])
    model = pa.make_model()
    assert(len(model.rules)==3)
    assert(len(model.monomers)==2)

def test_modification_assembly():
    classes = AddModification.__subclasses__() + \
              RemoveModification.__subclasses__()
    for mod_class in classes:
        _check_mod_assembly(mod_class)

def test_rule_annotation():
    a = Agent('A', db_refs={'HGNC': '1234'})
    b = Agent('B', db_refs={'HGNC': '5678'})

    def check_rule_annotation(stmt, policy):
        pa = PysbAssembler(policies=policy)
        pa.add_statements([stmt])
        model = pa.make_model()
        subj = [ann.object for ann in model.annotations
                if ann.predicate == 'rule_has_subject']
        obj = [ann.object for ann in model.annotations
                if ann.predicate == 'rule_has_object']
        assert len(subj) == 1
        assert subj[0] == 'A'
        assert len(obj) == 1
        assert obj[0] == 'B'

    classes = AddModification.__subclasses__() + \
              RemoveModification.__subclasses__()
    for mod_class in classes:
        stmt = mod_class(a, b)
        check_rule_annotation(stmt, 'one_step')
        check_rule_annotation(stmt, 'two_step')

    # Check ATP dependent phosphorylation
    stmt = Phosphorylation(a, b)
    check_rule_annotation(stmt, 'atp_dependent')
    stmt = Activation(a, b)
    check_rule_annotation(stmt, 'one_step')
    #Skip Autophosphorylation and Transphosphorylation for now
    #RasGef
    #RasGap

def test_activeform_site():
    a = Agent('A', db_refs={'HGNC': '1234'})
    b = Agent('B', db_refs={'HGNC': '5678'})
    b_phos = Agent('B', mods=[ModCondition('phosphorylation', 'Y', '200')],
                   db_refs={'HGNC': '5678'})
    st1 = Phosphorylation(a, b, 'S', '100')
    st2 = ActiveForm(b_phos, 'kinase', True)
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([st1, st2])
    model = pa.make_model()

# TODO Do the same for mutation condition
# TODO Localization condition
# TODO Bound condition
# TODO Unphosphorylated/unmodified forms (try ubiquitinated/acetylated lysine)

def test_activation_subj1():
    """No subject activity is defined."""
    st = Activation(Agent('a'), Agent('b'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(pa.model.monomers['a'].sites == [])
    left = pa.model.rules[0].reactant_pattern
    subj_left = left.complex_patterns[0].monomer_patterns[0]
    right = pa.model.rules[0].product_pattern
    subj_right = right.complex_patterns[0].monomer_patterns[0]
    assert(subj_left.site_conditions == {})
    assert(subj_right.site_conditions == {})

def test_activation_subj2():
    """Subject activity is defined explicitly."""
    a_act = Agent('a', activity=ActivityCondition('activity', True))
    st = Activation(a_act, Agent('b'))
    st2 = ActiveForm(Agent('a', mods=[ModCondition('phosphorylation')]),
                     'activity', True)
    pa = PysbAssembler()
    pa.add_statements([st, st2])
    pa.make_model()
    assert(pa.model.monomers['a'].sites == ['phospho'])
    left = pa.model.rules[0].reactant_pattern
    subj_left = left.complex_patterns[0].monomer_patterns[0]
    right = pa.model.rules[0].product_pattern
    subj_right = right.complex_patterns[0].monomer_patterns[0]
    assert(subj_left.site_conditions == {u'phospho': (u'p', WILD)})
    assert(subj_right.site_conditions == {u'phospho': (u'p', WILD)})

def test_activation_subj3():
    """Subject activity is defined implicitly by another statement."""
    a_act = Agent('a', activity=ActivityCondition('activity', True))
    st = Activation(a_act, Agent('b'))
    st2 = Activation(Agent('c'), Agent('a'))
    pa = PysbAssembler()
    pa.add_statements([st, st2])
    pa.make_model()
    assert(len(pa.model.rules) == 2)
    assert(pa.model.monomers['a'].sites == ['activity'])
    left = pa.model.rules[0].reactant_pattern
    subj_left = left.complex_patterns[0].monomer_patterns[0]
    right = pa.model.rules[0].product_pattern
    subj_right = right.complex_patterns[0].monomer_patterns[0]
    assert(subj_left.site_conditions == {u'activity': (u'active')})
    assert(subj_right.site_conditions == {u'activity': (u'active')})

def test_activation_subj4():
    """Subject activity is defined both explicitly and implicitly."""
    a_act = Agent('a', activity=ActivityCondition('activity', True))
    st = Activation(a_act, Agent('b'))
    st2 = Activation(Agent('c'), Agent('a'))
    st3 = ActiveForm(Agent('a', mods=[ModCondition('phosphorylation')]),
                     'activity', True)
    pa = PysbAssembler()
    pa.add_statements([st, st2, st3])
    pa.make_model()
    assert(set(pa.model.monomers['a'].sites) == set(['activity', 'phospho']))
    left = pa.model.rules[0].reactant_pattern
    subj_left = left.complex_patterns[0].monomer_patterns[0]
    right = pa.model.rules[0].product_pattern
    subj_right = right.complex_patterns[0].monomer_patterns[0]
    assert(subj_left.site_conditions == {u'phospho': (u'p', WILD)})
    assert(subj_right.site_conditions == {u'phospho': (u'p', WILD)})

def test_pysb_preassembler_replace_activities1():
    st1 = ActiveForm(Agent('a', location='nucleus'), 'activity', True)
    st2 = Phosphorylation(Agent('a', activity=ActivityCondition('activity', True)), Agent('b'))
    ppa = PysbPreassembler([st1, st2])
    ppa.replace_activities()
    assert(len(ppa.statements) == 2)
    assert(ppa.statements[1].enz.location == 'nucleus')

def test_pysb_preassembler_replace_activities2():
    a_act = Agent('a', activity=ActivityCondition('activity', True))
    st = Activation(a_act, Agent('b'))
    st2 = Activation(Agent('c'), Agent('a'))
    ppa = PysbPreassembler([st, st2])
    ppa.replace_activities()
    assert(len(ppa.statements) == 2)

def test_pysb_preassembler_replace_activities3():
    p = Agent('PPP2CA')
    bc = BoundCondition(p, False)
    erk = Agent('ERK')
    mek1 = Agent('MEK', mods=[ModCondition('phosphorylation',
                                           None, None, True)])
    mek2 = Agent('MEK', activity=ActivityCondition('activity', True),
                 bound_conditions=[bc])
    st2 = ActiveForm(mek1, 'activity', True)
    st1 = Phosphorylation(mek2, erk)
    ppa = PysbPreassembler([st1, st2])
    ppa.replace_activities()
    assert(len(ppa.statements) == 2)
    assert(ppa.statements[0].enz.mods)
    assert(ppa.statements[0].enz.bound_conditions)
