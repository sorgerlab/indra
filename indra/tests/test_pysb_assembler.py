from indra.pysb_assembler import PysbAssembler
from indra.statements import *

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
    print model.rules
    assert(len(model.rules)==2)
    assert(len(model.monomers)==3)

def test_pysb_assembler_phos1():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_phos2():
    hras = Agent('HRAS')
    enz = Agent('BRAF', bound_conditions=[BoundCondition(hras, True)])
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222')
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
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222')
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
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==4)

def test_pysb_assembler_autophos1():
    enz = Agent('MEK1')
    stmt = Autophosphorylation(enz, 'PhosphorylationSerine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_autophos2():
    raf1 = Agent('RAF1')
    enz = Agent('MEK1', bound_conditions=[BoundCondition(raf1, True)])
    stmt = Autophosphorylation(enz, 'PhosphorylationSerine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_autophos3():
    egfr = Agent('EGFR')
    enz = Agent('EGFR', bound_conditions=[BoundCondition(egfr, True)])
    stmt = Autophosphorylation(enz, 'PhosphorylationTyrosine', None)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_transphos1():
    egfr = Agent('EGFR')
    enz = Agent('EGFR', bound_conditions=[BoundCondition(egfr, True)])
    stmt = Transphosphorylation(enz, 'PhosphorylationTyrosine', None)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_actact1():
    egfr = Agent('EGFR')
    subj = Agent('GRB2', bound_conditions=[BoundCondition(egfr, True)])
    obj = Agent('SOS1')
    stmt = ActivityActivity(subj, 'act', 'increase', obj, 'act')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==3)

def test_pysb_assembler_dephos1():
    phos = Agent('PP2A')
    sub = Agent('MEK1')
    stmt = Dephosphorylation(phos, sub, 'PhosphorylationSerine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_dephos2():
    phos = Agent('PP2A')
    raf1 = Agent('RAF1')
    sub = Agent('MEK1', bound_conditions=[BoundCondition(raf1, True)])
    stmt = Dephosphorylation(phos, sub, 'PhosphorylationSerine', '222')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==3)

def test_pysb_assembler_rasgef1():
    gef = Agent('SOS1')
    ras = Agent('HRAS')
    stmt = RasGef(gef, 'catalytic', ras)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_rasgap1():
    gap = Agent('NF1')
    ras = Agent('HRAS')
    stmt = RasGap(gap, 'catalytic', ras)
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_actmod1():
    mek = Agent('MEK')
    erk = Agent('ERK')
    stmts = []
    stmts.append(ActivityModification(mek, ['PhosphorylationSerine', 
                                      'PhosphorylationSerine'], [218,222],
                                      'increases', 'act'))
    stmts.append(Phosphorylation(mek, erk, 'PhosphorylationThreonine', '185'))
    stmts.append(Phosphorylation(mek, erk, 'PhosphorylationTyrosine', '187'))
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==2)
    assert(len(model.monomers)==2)

def test_pysb_assembler_actmod2():
    mek = Agent('MEK')
    erk = Agent('ERK')
    stmts = []
    stmts.append(ActivityModification(mek, ['PhosphorylationSerine'], 
        [218], 'increases', 'act'))
    stmts.append(ActivityModification(mek, ['PhosphorylationSerine'], 
        [222], 'increases', 'act'))
    stmts.append(Phosphorylation(mek, erk, 'PhosphorylationThreonine', '185'))
    stmts.append(Phosphorylation(mek, erk, 'PhosphorylationTyrosine', '187'))
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==4)
    assert(len(model.monomers)==2)

def test_pysb_assembler_phos_twostep1():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222')
    pa = PysbAssembler(policies='two_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==3)
    assert(len(model.monomers)==2)

def test_pysb_assembler_dephos_twostep1():
    phos = Agent('PP2A')
    sub = Agent('MEK1')
    stmt = Dephosphorylation(phos, sub, 'PhosphorylationSerine', '222')
    pa = PysbAssembler(policies='two_step')
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==3)
    assert(len(model.monomers)==2)

def test_statement_specific_policies():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    phos = Agent('PP2A')
    stmt1 = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222')
    stmt2 = Dephosphorylation(phos, sub, 'PhosphorylationSerine', '222')
    policies = {'Phosphorylation': 'two_step',
                'Dephosphorylation': 'interactions_only'}
    pa = PysbAssembler(policies=policies)
    pa.add_statements([stmt1, stmt2])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==4)
    assert(len(model.monomers)==3)

def test_unspecified_statement_policies():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    phos = Agent('PP2A')
    stmt1 = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222')
    stmt2 = Dephosphorylation(phos, sub, 'PhosphorylationSerine', '222')
    policies = {'Phosphorylation': 'two_step',
                'other': 'interactions_only'}
    pa = PysbAssembler(policies=policies)
    pa.add_statements([stmt1, stmt2])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==4)
    assert(len(model.monomers)==3)

