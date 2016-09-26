from __future__ import absolute_import, print_function, unicode_literals
import os
from os.path import dirname, join
import sys
from indra import trips
import indra.statements
from indra.assemblers import PysbAssembler

def test_bind():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'The receptor tyrosine kinase EGFR binds the growth factor ligand EGF.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(has_hgnc_ref(st.members[0]))
    assert(has_hgnc_ref(st.members[1]))
    os.remove(fname)

def test_complex_bind():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'The EGFR-EGF complex binds another EGFR-EGF complex.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(has_hgnc_ref(st.members[0]))
    assert(has_hgnc_ref(st.members[1]))
    assert(st.members[0].bound_conditions)
    assert(st.members[1].bound_conditions)
    os.remove(fname)

def test_complex_bind2():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'The EGFR-EGFR complex binds GRB2.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'EGFR')
    assert(st.members[1].name == 'GRB2')
    assert(len(st.members[0].bound_conditions) == 1)
    assert(st.members[0].bound_conditions[0].agent.name == 'EGFR')
    os.remove(fname)

def test_complex_bind3():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'RAF binds to the RAS-GTP complex.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'RAF')
    assert(st.members[1].name == 'RAS')
    assert(len(st.members[1].bound_conditions) == 1)
    assert(st.members[1].bound_conditions[0].agent.name == 'GTP')
    os.remove(fname)

def test_complex_bind4():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'The RAF-RAS complex binds another RAF-RAS complex.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'RAF')
    assert(st.members[1].name == 'RAF')
    assert(len(st.members[0].bound_conditions) == 1)
    assert(st.members[0].bound_conditions[0].agent.name == 'RAS')
    assert(len(st.members[1].bound_conditions) == 1)
    assert(st.members[1].bound_conditions[0].agent.name == 'RAS')
    os.remove(fname)

def test_bound_mod():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'The adaptor protein GRB2 can bind EGFR that is phosphorylated on tyrosine.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(has_hgnc_ref(st.members[0]))
    assert(has_hgnc_ref(st.members[1]))
    assert(st.members[1].mods)
    assert(st.members[1].mods[0].mod_type == 'phosphorylation')
    assert(st.members[1].mods[0].residue == 'Y')
    os.remove(fname)

def test_not_bound_to():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'BRAF that is not bound to Vemurafenib binds MEK1.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'BRAF')
    assert(has_hgnc_ref(st.members[0]))
    assert(st.members[1].name == 'MAP2K1')
    assert(has_hgnc_ref(st.members[1]))
    assert(len(st.members[0].bound_conditions) == 1)
    assert(st.members[0].bound_conditions[0].agent.name.lower() == 'vemurafenib')
    assert(st.members[0].bound_conditions[0].is_bound == False)
    os.remove(fname)

def test_not_bound_to2():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'BRAF, not bound to Vemurafenib, phosphorylates MEK1.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_phosphorylation(st))
    assert(st.enz is not None)
    assert(st.sub is not None)
    assert(st.enz.name == 'BRAF')
    assert(has_hgnc_ref(st.enz))
    assert(st.sub.name == 'MAP2K1')
    assert(has_hgnc_ref(st.sub))
    assert(len(st.enz.bound_conditions) == 1)
    assert(st.enz.bound_conditions[0].agent.name.lower() == 'vemurafenib')
    assert(st.enz.bound_conditions[0].is_bound == False)
    os.remove(fname)

def test_not_bound_to3():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'SOS1 bound to GRB2 binds NRAS that is not bound to BRAF and GTP.' 
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(st.members[0].name == 'SOS1')
    assert(st.members[1].name == 'NRAS')
    assert(len(st.members[0].bound_conditions) == 1)
    assert(len(st.members[1].bound_conditions) == 2)
    os.remove(fname)

def test_not_bound_to4():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'BRAF that is not bound to NRAS and Vemurafenib binds BRAF ' +\
          'that is not bound to Vemurafenib.' 
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(st.members[0].name == 'BRAF')
    assert(st.members[1].name == 'BRAF')
    assert(len(st.members[0].bound_conditions) == 2)
    assert(len(st.members[1].bound_conditions) == 1)
    os.remove(fname)

def test_bound_to():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'NRAS, bound to GTP, binds BRAF.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'NRAS')
    assert(has_hgnc_ref(st.members[0]))
    assert(st.members[1].name == 'BRAF')
    assert(has_hgnc_ref(st.members[1]))
    assert(len(st.members[0].bound_conditions) == 1)
    assert(st.members[0].bound_conditions[0].agent.name == 'GTP')
    os.remove(fname)

def test_bound_to2():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'EGFR-bound GRB2 binds SOS1.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'GRB2')
    assert(has_hgnc_ref(st.members[0]))
    assert(st.members[1].name == 'SOS1')
    assert(has_hgnc_ref(st.members[1]))
    assert(len(st.members[0].bound_conditions) == 1)
    assert(st.members[0].bound_conditions[0].agent.name == 'EGFR')
    os.remove(fname)

def test_bound_to3():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'SOS1, bound to GRB2 binds RAS.' 
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'SOS1')
    assert(has_hgnc_ref(st.members[0]))
    assert(st.members[1].name == 'RAS')
    assert(len(st.members[0].bound_conditions) == 1)
    assert(st.members[0].bound_conditions[0].agent.name == 'GRB2')
    os.remove(fname)

def test_bound_to4(): 
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'RAS, bound to SOS1, binds GTP.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'RAS')
    assert(st.members[1].name == 'GTP')
    assert(len(st.members[0].bound_conditions) == 1)
    assert(st.members[0].bound_conditions[0].agent.name == 'SOS1')
    os.remove(fname)

def test_bound_to5():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'BRAF that is bound to NRAS and Vemurafenib binds ' +\
          'BRAF that is bound to NRAS and Vemurafenib.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_complex(st))
    assert(len(st.members) == 2)
    assert(st.members[0].name == 'BRAF')
    assert(st.members[1].name == 'BRAF')
    assert(len(st.members[0].bound_conditions) == 2)
    assert(len(st.members[1].bound_conditions) == 2)
    os.remove(fname)

def test_transphosphorylate():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'EGFR, bound to EGFR, transphosphorylates itself.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_transphosphorylation(st))
    assert(st.enz is not None)
    assert(st.enz.name == 'EGFR')
    assert(has_hgnc_ref(st.enz))
    assert(st.residue == None)
    assert(len(st.enz.bound_conditions) == 1)
    assert(st.enz.bound_conditions[0].agent.name == 'EGFR')
    os.remove(fname)

def test_transphosphorylate2():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'EGFR, bound to EGFR, transphosphorylates itself on a tyrosine residue.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_transphosphorylation(st))
    assert(st.enz is not None)
    assert(st.enz.name == 'EGFR')
    assert(has_hgnc_ref(st.enz))
    assert(st.residue == 'Y')
    assert(len(st.enz.bound_conditions) == 1)
    assert(st.enz.bound_conditions[0].agent.name == 'EGFR')
    os.remove(fname)

def test_act_mod():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'MEK1, phosphorylated at Ser218 and Ser222, is activated.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_actmod(st))
    assert(st.agent is not None)
    assert(st.agent.name == 'MAP2K1')
    residues = [m.residue for m in st.agent.mods]
    positions = [m.position for m in st.agent.mods]
    assert(residues == ['S', 'S'])
    assert(positions == ['218', '222'])
    assert(st.is_active)
    os.remove(fname)

def test_bound_phosphorylate():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'RAF, bound to RAF, phosphorylates MEK1.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_phosphorylation(st))
    assert(st.enz is not None)
    assert(st.enz.name == 'RAF')
    assert(st.sub is not None)
    assert(st.sub.name == 'MAP2K1')
    assert(st.residue == None)
    os.remove(fname)

'''
def test_bound_not_bound_phosphorylate():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'BRAF-bound BRAF that is not bound to Vemurafenib phosphorylates MEK1.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_phosphorylation(st))
    assert(st.enz is not None)
    assert(st.enz.name == 'MAP2K1')
    assert(st.monomer.mod == ['PhosphorylationSerine', 'PhosphorylationSerine'])
    assert(st.monomer.mod_pos == ['218', '222'])
    assert(st.relationship == 'increases')
    os.remove(fname)
'''

def test_act_phosphorylate():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'Active MEK1 phosphorylates ERK2 at Tyr187.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_phosphorylation(st))
    assert(st.enz is not None)
    assert(st.enz.name == 'MAP2K1')
    assert(st.sub is not None)
    assert(st.sub.name == 'MAPK1')
    assert(st.residue == 'Y')
    assert(st.position == '187')
    os.remove(fname)

def test_dephosphorylate():
    fname = sys._getframe().f_code.co_name + '.xml'
    txt = 'DUSP6 dephosphorylates ERK2 at Tyr187.'
    tp = trips.process_text(txt, fname, False)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(is_dephosphorylation(st))
    assert(st.enz is not None)
    assert(st.enz.name == 'DUSP6')
    assert(st.sub is not None)
    assert(st.sub.name == 'MAPK1')
    assert(st.residue == 'Y')
    assert(st.position == '187')
    os.remove(fname)

def is_complex(statement):
    return isinstance(statement, indra.statements.Complex)

def is_phosphorylation(statement):
    return isinstance(statement, indra.statements.Phosphorylation)

def is_transphosphorylation(statement):
    return isinstance(statement, indra.statements.Transphosphorylation)

def is_actmod(statement):
    return isinstance(statement, indra.statements.ActiveForm)

def is_dephosphorylation(statement):
    return isinstance(statement, indra.statements.Dephosphorylation)

def has_hgnc_ref(agent):
    return ('HGNC' in agent.db_refs.keys())
