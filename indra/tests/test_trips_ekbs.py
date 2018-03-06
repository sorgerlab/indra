from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import os
from indra.statements import *
from indra.sources import trips

path_this = os.path.dirname(os.path.abspath(__file__))

def assert_if_hgnc_then_up(st):
    agents = st.agent_list()
    for a in agents:
        if a is not None:
            up_id = a.db_refs.get('UP')
            hgnc_id = a.db_refs.get('HGNC')
            if hgnc_id and not up_id:
                assert(False)

def assert_grounding_value_or_none(st):
    agents = st.agent_list()
    for a in agents:
        if a is not None:
            for k, v in a.db_refs.items():
                # Make sure there are no empty strings/lists
                if not v:
                    assert(v is None)


def process_sentence_xml(sentence):
    fname = re.sub('[^a-zA-Z0-9]', '_', sentence[:-1]) + '.ekb'
    path = os.path.join(path_this, 'trips_ekbs', fname)
    with open(path, 'rb') as fh:
        xml = fh.read().decode('utf-8')
    tp = trips.process_xml(xml)
    return tp

def assert_onestmt(tp):
    assert(tp is not None)
    assert(len(tp.statements) == 1)

def assert_evidence(stmt):
    assert(len(stmt.evidence) == 1)
    assert(stmt.evidence[0].source_api == 'trips')
    assert(stmt.evidence[0].text)

def assert_modtype(stmt, mod_type):
    assert(isinstance(stmt, mod_type))
    assert(stmt.enz is not None)
    assert(stmt.sub is not None)

def test_1():
    sentence = 'The receptor tyrosine kinase EGFR binds the growth factor ligand EGF.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Complex))
    assert(len(st.members) == 2)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_2():
    sentence = 'The EGFR-EGF complex binds another EGFR-EGF complex.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Complex))
    assert(len(st.members) == 2)
    e1, e2 = st.members
    assert(e1.bound_conditions[0].is_bound)
    assert(e2.bound_conditions[0].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_3():
    sentence = 'The EGFR-EGFR complex binds GRB2.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Complex))
    assert(len(st.members) == 2)
    e, g = st.members
    assert(e.bound_conditions[0].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_4():
    sentence = 'EGFR-bound GRB2 binds SOS1.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Complex))
    assert(len(st.members) == 2)
    g, s = st.members
    assert(g.bound_conditions[0].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_5():
    sentence = 'GRB2-bound SOS1 binds NRAS that is not bound to BRAF.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Complex))
    assert(len(st.members) == 2)
    s, n = st.members
    assert(s.bound_conditions[0].is_bound)
    assert(not n.bound_conditions[0].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_6():
    sentence = 'SOS1-bound NRAS binds GTP.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Complex))
    assert(len(st.members) == 2)
    n, g = st.members
    assert(n.bound_conditions[0].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_7():
    sentence = 'GTP-bound NRAS that is not bound to SOS1 binds BRAF.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Complex))
    assert(len(st.members) == 2)
    n, b = st.members
    assert(len(n.bound_conditions) == 2)
    assert(n.bound_conditions[0].is_bound)
    assert(not n.bound_conditions[1].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_8():
    sentence = 'Vemurafenib binds BRAF.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Complex))
    assert(len(st.members) == 2)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_9():
    sentence = 'BRAF V600E that is not bound to Vemurafenib phosphorylates MAP2K1.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Phosphorylation))
    assert(st.enz is not None)
    assert(st.sub is not None)
    assert(st.enz.mutations)
    assert(not st.enz.bound_conditions[0].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_10():
    sentence = 'PP2A-alpha dephosphorylates MAP2K1 that is not bound to ERK2.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Dephosphorylation))
    assert(st.enz is not None)
    assert(st.sub is not None)
    assert(not st.sub.bound_conditions[0].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_11():
    sentence = 'Phosphorylated MAP2K1 is activated.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ActiveForm))
    assert(st.agent is not None)
    assert(st.agent.mods[0].mod_type == 'phosphorylation')
    assert(st.agent.mods[0].is_modified)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_12():
    sentence = 'Active MAP2K1 that is not bound to PP2A-alpha phosphorylates ERK2.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Phosphorylation))
    assert(st.enz is not None)
    assert(st.sub is not None)
    assert(st.enz.activity.activity_type == 'activity')
    assert(not st.enz.bound_conditions[0].is_bound)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_13():
    sentence = 'DUSP6 dephosphorylates ERK2.'
    tp = process_sentence_xml(sentence)
    assert(tp is not None)
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, Dephosphorylation))
    assert(st.enz is not None)
    assert(st.sub is not None)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert_evidence(st)

def test_14():
    sentence = 'MAP2K1 ubiquitinates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Ubiquitination)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_15():
    sentence = 'MAP2K1 ribosylates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Ribosylation)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)

def test_16():
    sentence = 'MAP2K1 hydroxylates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Hydroxylation)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)

def test_17():
    sentence = 'MAP2K1 acetylates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Acetylation)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_18():
    sentence = 'MAP2K1 farnesylates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Farnesylation)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_19():
    sentence = 'MAP2K1 deubiquitinates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Deubiquitination)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_20():
    sentence = 'MAP2K1 deribosylates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Deribosylation)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_21():
    sentence = 'MAP2K1 dehydroxylates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Dehydroxylation)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_22():
    sentence = 'MAP2K1 deacetylates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Deacetylation)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)

def test_23():
    sentence = 'MAP2K1 defarnesylates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert_modtype(st, Defarnesylation)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_24():
    sentence = 'Ubiquitinated MAPK1 is degraded.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, DecreaseAmount))
    assert(st.subj is None)
    assert(st.obj is not None)
    assert(len(st.obj.mods) == 1)
    assert(st.obj.mods[0].mod_type == 'ubiquitination')
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_25():
    sentence = 'MAPK1 is synthesized.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, IncreaseAmount))
    assert(st.subj is None)
    assert(st.obj is not None)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_26():
    sentence = 'MAP2K1 transcribes MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, IncreaseAmount))
    assert(st.subj is not None)
    assert(st.obj is not None)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_27():
    sentence = 'MAP2K1 synthesizes MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, IncreaseAmount))
    assert(st.subj is not None)
    assert(st.obj is not None)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_28():
    sentence = 'MAP2K1 degrades MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, DecreaseAmount))
    assert(st.subj is not None)
    assert(st.obj is not None)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_29():
    sentence = 'MAPK1 translocates to the nucleus.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Translocation))
    assert(st.agent is not None)
    assert(st.to_location == 'nucleus')
    assert(st.from_location is None)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_30():
    sentence = 'MAPK1 translocates from the nucleus.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Translocation))
    assert(st.agent is not None)
    assert(st.from_location == 'nucleus')
    assert(st.to_location is None)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_31():
    sentence = 'MAPK1 translocates from the plasma membrane to the nucleus.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Translocation))
    assert(st.agent is not None)
    assert(st.to_location == 'nucleus')
    assert(st.from_location == 'plasma membrane')
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_32():
    sentence = 'EGF leads to the activation of MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Activation))
    assert(st.obj is not None)
    assert(st.subj is not None)
    assert(st.is_activation)
    assert(not st.evidence[0].epistemics['direct'])
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

'''
TODO: put back when this is implemented
def test_33():
    sentence = 'Stimulation by EGF activates MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Activation))
    assert(st.obj is not None)
    assert(st.subj is not None)
    assert(st.is_activation)
    assert(not st.evidence[0].epistemics['direct'])
    assert_evidence(st)
'''

def test_34():
    sentence = 'Vemurafenib leads to the deactivation of MAPK1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Inhibition))
    assert(st.obj is not None)
    assert(st.subj is not None)
    assert(not st.is_activation)
    assert(not st.evidence[0].epistemics['direct'])
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_35():
    sentence = 'EGFR autophosphorylates itself.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Autophosphorylation))
    assert(st.enz is not None)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_36():
    sentence = 'EGFR autophosphorylates itself on Y1234.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Autophosphorylation))
    assert(st.enz is not None)
    assert(st.residue == 'Y')
    assert(st.position == '1234')
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_37():
    sentence = 'EGFR bound to EGFR transphosphorylates itself.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Transphosphorylation))
    assert(st.enz is not None)
    assert(st.enz.bound_conditions)
    assert_evidence(st)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)

def test_38():
    sentence = 'TCRA activates NEDD4, MEK1, CK2, PIP3 and mTORC2.'
    tp = process_sentence_xml(sentence)
    assert(len(tp.statements) == 5)
    for st in tp.statements:
        assert(isinstance(st, Activation))
        assert_grounding_value_or_none(st)
        # Here we don't assert if_hgnc_then_up
        # because TCRA maps to an HGNC gene that has
        # no corresponding protein
        #assert_if_hgnc_then_up(st)

def test_39():
    sentence = 'FGF2 activates PI3K/Akt/mTOR and MAPK/ERK.'
    tp = process_sentence_xml(sentence)
    # For now, this should not return any statements
    assert(not tp.statements)

def test_40():
    sentence = 'Ras-activated SAF-1 that binds to a bona fide SAF-1-binding element.'
    tp = process_sentence_xml(sentence)
    # All events here are static so nothing should be extracted
    assert(not tp.statements)

def test_41():
    sentence = 'NFKB synthesizes IKB in the nucleus.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(st.obj.location is not None)

def test_42():
    sentence = 'RAF1 activates MAP2K1.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    raf1 = st.subj
    assert(raf1.name == 'RAF1')
    assert(raf1.db_refs.get('HGNC') == '9829')

def test_43():
    sentence = 'Phosphorylated ERK is active.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    erk = st.agent
    assert(erk.name == 'ERK')
    assert(erk.db_refs.get('FPLX') == 'ERK')
    assert(len(erk.mods) == 1)

def test_44():
    sentence = 'p53 positively regulates the transcription of mdm2.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, IncreaseAmount))
    p53 = st.subj
    mdm2 = st.obj
    assert(p53.name == 'TP53')
    assert(mdm2.name == 'MDM2')

'''
# Not sure if this is good to extract like this so leaving out for now
def test_45():
    sentence = 'p53 increases mdm2.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, IncreaseAmount))
    p53 = st.subj
    mdm2 = st.obj
    assert(p53.name == 'TP53')
    assert(mdm2.name == 'MDM2')
'''

def test_46():
    sentence = 'p53 increases the transcription of mdm2.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, IncreaseAmount))
    p53 = st.subj
    mdm2 = st.obj
    assert(p53.name == 'TP53')
    assert(mdm2.name == 'MDM2')

def test_47():
    sentence = 'p53 decreases the transcription of mdm2.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, DecreaseAmount))
    p53 = st.subj
    mdm2 = st.obj
    assert(p53.name == 'TP53')
    assert(mdm2.name == 'MDM2')

def test_48():
    sentence = 'p53 downregulates the transcription of mdm2.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, DecreaseAmount))
    p53 = st.subj
    mdm2 = st.obj
    assert(p53.name == 'TP53')
    assert(mdm2.name == 'MDM2')

def test_49():
    sentence = 'Ras converts GTP into GDP.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Conversion))
    ras = st.subj
    gtp = st.obj_from[0]
    gdp = st.obj_to[0]
    assert(ras.name == 'RAS')
    assert(gtp.name == 'GTP')
    assert(gdp.name == 'GDP')

def test_50():
    sentence = 'GTP is converted into GDP.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Conversion))
    ras = st.subj
    gtp = st.obj_from[0]
    gdp = st.obj_to[0]
    assert(ras is None)
    assert(gtp.name == 'GTP')
    assert(gdp.name == 'GDP')

def test_51():
    sentence = 'RAS converts GTP into GDP and GMP.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Conversion))
    ras = st.subj
    gtp = st.obj_from[0]
    gdp_gmp = st.obj_to
    assert(ras.name == 'RAS')
    assert(gtp.name == 'GTP')
    assert(gdp_gmp[0].name == 'GDP')
    assert(gdp_gmp[1].name == 'GMP')

def test_52():
    sentence = 'PTEN catalyzes the conversion of PIP3 to PIP2.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Conversion))
    pten = st.subj
    pip3 = st.obj_from[0]
    pip2 = st.obj_to[0]
    assert(pten.name == 'PTEN')
    assert(pip2.name.startswith('PIP'))
    assert(pip2.name.endswith('2'))
    assert(pip3.name.startswith('PIP'))
    assert(pip3.name.endswith('3'))

def test_53():
    sentence = 'MEK increases the phosphorylation of ERK.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Phosphorylation))
    mek = st.enz
    erk = st.sub
    assert(mek.name == 'MEK')
    assert(erk.name == 'ERK')
    for ev in st.evidence:
        assert ev.epistemics.get('direct') == False

def test_54():
    sentence = 'EGF leads to the phosphorylation of ERK.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Phosphorylation))
    mek = st.enz
    erk = st.sub
    assert(mek.name == 'EGF')
    assert(erk.name == 'ERK')
    for ev in st.evidence:
        assert ev.epistemics.get('direct') == False

def test_55():
    sentence = 'Unphosphorylated ERK is degraded.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, DecreaseAmount))
    erk = st.obj
    assert(erk.name == 'ERK')
    assert len(erk.mods) == 1
    assert erk.mods[0].is_modified == False

def test_56():
    sentence = 'Activated TGFBR1 phosphorylates SMURF2.'
    tp = process_sentence_xml(sentence)
    assert_onestmt(tp)
    st = tp.statements[0]
    assert(isinstance(st, Phosphorylation))
    tgfbr1 = st.enz
    assert(tgfbr1.name == 'TGFBR1')
    assert tgfbr1.activity is not None
    assert tgfbr1.activity.activity_type == 'activity'
    assert tgfbr1.activity.is_active == True
