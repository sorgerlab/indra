from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers import PysbAssembler
from indra.assemblers.pysb_assembler import get_agent_rule_str, _n
from indra.statements import *
from pysb import bng

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
    stmt = Activation(subj, 'activity', obj, 'activity', True)
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
    stmt = RasGef(gef, 'catalytic', ras)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_rasgap1():
    gap = Agent('NF1')
    ras = Agent('HRAS')
    stmt = RasGap(gap, 'catalytic', ras)
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

def test_pysb_assembler_actmod2():
    mek = Agent('MEK')
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
    stmt = Activation(subj, 'activity', obj, 'activity', True)
    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_activity_activity():
    subj = Agent('KRAS')
    obj = Agent('BRAF')
    stmt = Activation(subj, 'activity', obj, 'activity', True)
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_activity_activity2():
    subj = Agent('Vemurafenib')
    obj = Agent('BRAF')
    stmt = Activation(subj, None, obj, 'activity', False)
    pa = PysbAssembler(policies='interactions_only')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_activity_activity3():
    subj = Agent('Vemurafenib')
    obj = Agent('BRAF')
    stmt = Activation(subj, None, obj, 'activity', False)
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_rule_name_str_1():
    s = get_agent_rule_str(Agent('BRAF'))
    assert(s == 'BRAF')

def test_rule_name_str_2():
    a = Agent('GRB2',
              bound_conditions=[BoundCondition(Agent('EGFR'), True)])
    s = get_agent_rule_str(a)
    assert(s == 'GRB2_EGFR')

def test_rule_name_str_3():
    a = Agent('GRB2',
              bound_conditions=[BoundCondition(Agent('EGFR'), False)])
    s = get_agent_rule_str(a)
    assert(s == 'GRB2_nEGFR')

def test_rule_name_str_4():
    a = Agent('BRAF', mods=[ModCondition('phosphorylation', 'serine')])
    s = get_agent_rule_str(a)
    assert(s == 'BRAF_phosphoS')

def test_rule_name_str_5():
    a = Agent('BRAF', mods=[ModCondition('phosphorylation', 'serine', '123')])
    s = get_agent_rule_str(a)
    assert(s == 'BRAF_phosphoS123')

def test_neg_act_mod():
    mc = ModCondition('phosphorylation', 'serine', '123', False)
    st1 = ActiveForm(Agent('BRAF', mods=[mc]), 'active', True)
    st2 = Phosphorylation(Agent('BRAF'), Agent('MAP2K2'))
    pa = PysbAssembler(policies='one_step')
    pa.add_statements([st1, st2])
    pa.make_model()
    assert(len(pa.model.rules) == 1)
    r = pa.model.rules[0]
    braf = r.reactant_pattern.complex_patterns[0].monomer_patterns[0]
    assert(braf.monomer.name == 'BRAF')
    assert(braf.site_conditions == {'S123': 'u'})

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
    assert(braf.site_conditions == {'S123': 'p'})

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
    assert(braf.site_conditions == {'S123': 'u'})

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
    assert(pa.model.parameters['MAP2K1_0'].value < 1000)
    assert(pa.model.parameters['MAPK3_0'].value < 1000)
    pa.set_context('A375_SKIN')
    assert(pa.model.parameters['MAP2K1_0'].value > 10000)
    assert(pa.model.parameters['MAPK3_0'].value > 10000)

def test_set_context_monomer_notfound():
    st = Phosphorylation(Agent('MAP2K1'), Agent('XYZ'))
    pa = PysbAssembler()
    pa.add_statements([st])
    pa.make_model()
    assert(pa.model.parameters['MAP2K1_0'].value < 1000)
    assert(pa.model.parameters['XYZ_0'].value < 1000)
    pa.set_context('A375_SKIN')
    assert(pa.model.parameters['MAP2K1_0'].value > 10000)
    assert(pa.model.parameters['XYZ_0'].value < 1000)

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
    assert(len(pa.model.annotations) == 2)

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

def test_name_standardize():
    n = _n('.*/- ^&#@$')
    assert(isinstance(n, str))
    assert(n == '__________')
    n = _n('14-3-3')
    assert(isinstance(n, str))
    assert(n == 'p14_3_3')
    n = _n('\U0001F4A9bar')
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
