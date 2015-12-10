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

def test_pysb_assembler_complex2():
    member1 = Agent('BRAF', bound_to='HRAS')
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
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_phos2():
    enz = Agent('BRAF', bound_to='HRAS')
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==3)

def test_pysb_assembler_phos3():
    enz = Agent('BRAF', bound_to='HRAS')
    sub = Agent('MEK1', bound_to='ERK1')
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==4)

def test_pysb_assembler_phos4():
    enz = Agent('BRAF', bound_to='HRAS')
    sub = Agent('MEK1', bound_to='ERK1', bound_neg=True)
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==4)

def test_pysb_assembler_autophos1():
    enz = Agent('MEK1')
    stmt = Autophosphorylation(enz, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_autophos2():
    enz = Agent('MEK1', bound_to='RAF1')
    stmt = Autophosphorylation(enz, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_autophos3():
    enz = Agent('EGFR', bound_to='EGFR')
    stmt = Autophosphorylation(enz, 'PhosphorylationTyrosine', None, '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_transphos1():
    enz = Agent('EGFR', bound_to='EGFR')
    stmt = Transphosphorylation(enz, 'PhosphorylationTyrosine', None, '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==1)

def test_pysb_assembler_actact1():
    subj = Agent('GRB2', bound_to='EGFR')
    obj = Agent('SOS1')
    stmt = ActivityActivity(subj, 'act', 'increase', obj, 'act', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==3)

def test_pysb_assembler_dephos1():
    phos = Agent('PP2A')
    sub = Agent('MEK1')
    stmt = Dephosphorylation(phos, sub, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_dephos2():
    phos = Agent('PP2A')
    sub = Agent('MEK1', bound_to='RAF1')
    stmt = Dephosphorylation(phos, sub, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==3)

def test_pysb_assembler_rasgef1():
    gef = Agent('SOS1')
    ras = Agent('HRAS')
    stmt = RasGef(gef, 'catalytic', ras, '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==1)
    assert(len(model.monomers)==2)

def test_pysb_assembler_rasgap1():
    gap = Agent('NF1')
    ras = Agent('HRAS')
    stmt = RasGap(gap, 'catalytic', ras, '', '', '', '')
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
        'PhosphorylationSerine'], [218,222], 'increases', 'act', '', '', '', ''))
    stmts.append(Phosphorylation(mek, erk, 'PhosphorylationThreonine', '185', '', '', '', ''))
    stmts.append(Phosphorylation(mek, erk, 'PhosphorylationTyrosine', '187', '', '', '', ''))
    
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
        [218], 'increases', 'act', '', '', '', ''))
    stmts.append(ActivityModification(mek, ['PhosphorylationSerine'], 
        [222], 'increases', 'act', '', '', '', ''))
    stmts.append(Phosphorylation(mek, erk, 'PhosphorylationThreonine', '185', '', '', '', ''))
    stmts.append(Phosphorylation(mek, erk, 'PhosphorylationTyrosine', '187', '', '', '', ''))
    
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    print model.rules
    assert(len(model.rules)==4)
    assert(len(model.monomers)==2)

def test_pysb_assembler_phos_twostep1():
    enz = Agent('BRAF')
    sub = Agent('MEK1')
    stmt = Phosphorylation(enz, sub, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model(policies='two_step')
    print model.rules
    assert(len(model.rules)==3)
    assert(len(model.monomers)==2)

def test_pysb_assembler_dephos_twostep1():
    phos = Agent('PP2A')
    sub = Agent('MEK1')
    stmt = Dephosphorylation(phos, sub, 'PhosphorylationSerine', '222', '', '', '', '')
    pa = PysbAssembler()
    pa.add_statements([stmt])
    model = pa.make_model(policies='two_step')
    print model.rules
    assert(len(model.rules)==3)
    assert(len(model.monomers)==2)
