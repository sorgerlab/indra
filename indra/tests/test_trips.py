from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
from os.path import dirname, join
import indra.statements as ist
from indra import trips
from indra.assemblers import PysbAssembler
from indra.util import unicode_strs

test_small_file = join(dirname(__file__), 'test_small.xml')

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

def test_phosphorylation():
    tp = trips.process_text('BRAF phosphorylates MEK1 at Ser222.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.Phosphorylation))
    assert(st.residue == 'S')
    assert(st.position == '222')
    assert(st.evidence)
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert unicode_strs((tp, st))

def test_mod_cond():
    tp = trips.process_text('Phosphorylated BRAF binds ubiquitinated MAP2K1.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.Complex))
    braf = st.members[0]
    mek = st.members[1]
    assert(len(braf.mods) == 1)
    assert(braf.mods[0].mod_type == 'phosphorylation')
    assert(len(mek.mods) == 1)
    assert(mek.mods[0].mod_type == 'ubiquitination')
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_ubiquitination():
    tp = trips.process_text('MDM2 ubiquitinates TP53.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.Ubiquitination))
    assert unicode_strs((tp, st))
    assert_grounding_value_or_none(st)
    assert_if_hgnc_then_up(st)
    assert(st.evidence)

def test_phosphorylation_noresidue():
    tp = trips.process_text('BRAF phosphorylates MEK1.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.Phosphorylation))
    assert(st.residue is None)
    assert(st.position is None)
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_phosphorylation_nosite():
    tp = trips.process_text('BRAF phosphorylates MEK1 at Serine.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.Phosphorylation))
    assert(st.residue == 'S')
    assert(st.position is None)
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_actmod():
    tp = trips.process_text('MEK1 phosphorylated at Ser222 is activated.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.ActiveForm))
    assert(isinstance(st.agent.mods[0], ist.ModCondition))
    assert(st.agent.mods[0].mod_type == 'phosphorylation')
    assert(st.agent.mods[0].residue == 'S')
    assert(st.agent.mods[0].position == '222')
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_actmods():
    tp = trips.process_text('MEK1 phosphorylated at Ser 218 and Ser222 is activated.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.ActiveForm))
    assert(isinstance(st.agent.mods[0], ist.ModCondition))
    assert(isinstance(st.agent.mods[1], ist.ModCondition))
    assert(st.agent.mods[0].mod_type == 'phosphorylation')
    assert(st.agent.mods[0].residue == 'S')
    assert(st.agent.mods[0].position == '218')
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_actform_bound():
    tp = trips.process_text('HRAS bound to GTP is activated.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.ActiveForm))
    assert(isinstance(st.agent.bound_conditions[0], ist.BoundCondition))
    assert(st.agent.bound_conditions[0].agent.name == 'GTP')
    assert(st.agent.bound_conditions[0].is_bound == True)
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_actform_muts():
    tp = trips.process_text('BRAF V600E is activated.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.ActiveForm))
    assert(isinstance(st.agent.mutations[0], ist.MutCondition))
    assert(st.agent.mutations[0].residue_from == 'V')
    assert(st.agent.mutations[0].residue_to == 'E')
    assert(st.agent.mutations[0].position == '600')
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_actmods2():
    tp = trips.process_text('BRAF phosphorylated at Ser536 binds MEK1.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.Complex))
    braf = st.members[0]
    assert(braf.mods[0].mod_type == 'phosphorylation')
    assert(braf.mods[0].residue == 'S')
    assert(braf.mods[0].position == '536')
    assert unicode_strs((tp, st, braf))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_synthesis():
    tp = trips.process_text('NFKB transcribes IKB.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.IncreaseAmount))
    assert(st.subj is not None)
    assert(st.obj is not None)
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)

def test_degradation():
    tp = trips.process_text('MDM2 degrades TP53.')
    assert(len(tp.statements) == 1)
    st = tp.statements[0]
    assert(isinstance(st, ist.DecreaseAmount))
    assert(st.subj is not None)
    assert(st.obj is not None)
    assert unicode_strs((tp, st))
    assert_if_hgnc_then_up(st)
    assert_grounding_value_or_none(st)
    assert(st.evidence)
