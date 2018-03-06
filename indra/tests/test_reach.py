from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import unittest
from nose.plugins.attrib import attr
from indra.sources import reach
from indra.sources.reach.processor import ReachProcessor
from indra.util import unicode_strs
from indra.statements import IncreaseAmount, DecreaseAmount, Dephosphorylation

# Change this list to control what modes of
# reading are enabled in tests
offline_modes = [True]

def test_parse_site_text():
    text = ['threonine 185', 'thr 185', 'thr-185',
            'threonine residue 185', 'T185']
    assert unicode_strs(text)
    for t in text:
        residue, site = ReachProcessor._parse_site_text(t)
        assert(residue == 'T')
        assert(site == '185')
        assert unicode_strs((residue, site))

def test_parse_site_text_number():
    t = '135'
    residue, site = ReachProcessor._parse_site_text(t)
    assert(residue is None)
    assert(site == '135')
    assert(unicode_strs(site))

def test_parse_site_text_number_first():
    t = '293T'
    residue, site = ReachProcessor._parse_site_text(t)
    assert(residue == 'T')
    assert(site == '293')
    assert(unicode_strs((residue, site)))

def test_parse_site_text_number_first_space():
    t = '293 T'
    residue, site = ReachProcessor._parse_site_text(t)
    assert(residue == 'T')
    assert(site == '293')
    assert(unicode_strs((residue, site)))

def test_parse_site_text_other_aa():
    t = 'A431'
    residue, site = ReachProcessor._parse_site_text(t)
    assert(residue == 'A')
    assert(site == '431')
    assert(unicode_strs((residue, site)))

def test_parse_site_residue_only():
    text = ['serine residue', 'serine', 'a serine site', 's', 'ser']
    assert unicode_strs(text)
    for t in text:
        residue, site = ReachProcessor._parse_site_text(t)
        assert unicode_strs((residue, site))
        assert(residue == 'S')
        assert(site is None)

def test_valid_name():
    assert(ReachProcessor._get_valid_name('') == '')
    assert(ReachProcessor._get_valid_name('a') == 'a')
    assert(ReachProcessor._get_valid_name('Name123') == 'Name123')
    assert(ReachProcessor._get_valid_name('<>#~!,./][;-') == '____________')
    assert(ReachProcessor._get_valid_name('PI3 Kinase') == 'PI3_Kinase')
    assert(ReachProcessor._get_valid_name('14-3-3') == 'p14_3_3')

def test_phosphorylate():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 phosphorylates ERK2.', offline=offline)
        assert(len(rp.statements) == 1)
        s = rp.statements[0]
        assert (s.enz.name == 'MAP2K1')
        assert (s.sub.name == 'MAPK1')
        assert unicode_strs(rp.statements)


def test_indirect_phosphorylate():
    txt = 'DUSP decreases the phosphorylation of ERK.'
    for offline in offline_modes:
        rp = reach.process_text(txt, offline=offline)
        assert len(rp.statements) == 1
        s = rp.statements[0]
        assert isinstance(s, Dephosphorylation)
        assert s.enz.name == 'DUSP'
        assert s.sub.name == 'ERK'
        assert s.evidence[0].epistemics.get('direct') == False


def test_regulate_amount():
    for offline in offline_modes:
        rp = reach.process_text('ERK increases the transcription of DUSP.',
                                offline=offline)
        assert(len(rp.statements) == 1)
        s = rp.statements[0]
        assert(isinstance(s, IncreaseAmount))
        assert (s.subj.name == 'ERK')
        assert (s.obj.name == 'DUSP')
        assert unicode_strs(rp.statements)
        rp = reach.process_text('ERK decreases the amount of DUSP.',
                                offline=offline)
        assert(len(rp.statements) == 1)
        s = rp.statements[0]
        assert(isinstance(s, DecreaseAmount))
        assert (s.subj.name == 'ERK')
        assert (s.obj.name == 'DUSP')
        assert unicode_strs(rp.statements)

def test_multiple_enzymes():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 and MEK2 phosphorylate ERK1.',
                                offline=offline)
        assert(len(rp.statements) == 2)
        s = rp.statements[0]
        if s.enz.name == 'MAP2K1':
            assert(rp.statements[1].enz.name == 'MAP2K2')
        else:
            assert(rp.statements[1].enz.name == 'MAP2K1')
        assert (s.sub.name == 'MAPK3')
        s = rp.statements[1]
        assert (s.sub.name == 'MAPK3')
        assert unicode_strs(rp.statements)

def test_activate():
    for offline in offline_modes:
        rp = reach.process_text('HRAS activates BRAF.', offline=offline)
        assert(len(rp.statements) == 1)
        s = rp.statements[0]
        assert (s.subj.name == 'HRAS')
        assert (s.obj.name == 'BRAF')
        assert unicode_strs(rp.statements)

def test_bind():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 binds ERK2.', offline=offline)
        assert(len(rp.statements) == 1)
        assert unicode_strs(rp.statements)

def test_be_grounding():
    for offline in offline_modes:
        rp = reach.process_text('MEK activates ERK.', offline=offline)
        assert(len(rp.statements) == 1)
        assert unicode_strs(rp.statements)
        if offline == True:
            st = rp.statements[0]
            assert(st.subj.db_refs.get('FPLX') == 'MEK')
            assert(st.obj.db_refs.get('FPLX') == 'ERK')

def test_activity():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 activates ERK2.', offline=offline)
        assert(len(rp.statements) == 1)
        assert unicode_strs(rp.statements)

def test_mutation():
    for offline in offline_modes:
        rp = reach.process_text('BRAF(V600E) phosphorylates MEK.',
                                offline=offline)
        assert(len(rp.statements) == 1)
        braf = rp.statements[0].enz
        assert(braf.name == 'BRAF')
        assert(len(braf.mutations) == 1)
        assert(braf.mutations[0].position == '600')
        assert(braf.mutations[0].residue_from == 'V')
        assert(braf.mutations[0].residue_to == 'E')
        assert unicode_strs(rp.statements)

def test_process_unicode():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 binds ERK2\U0001F4A9.', offline=offline)
        assert unicode_strs(rp.statements)

@attr('slow')
def test_process_pmc():
    for offline in offline_modes:
        rp = reach.process_pmc('PMC4338247', offline=offline)
        for stmt in rp.statements:
            assert_pmid(stmt)
        assert unicode_strs(rp.statements)

def test_process_unicode_abstract():
    for offline in offline_modes:
        rp = reach.process_pubmed_abstract('27749056', offline=offline)
        assert unicode_strs(rp.statements)

def test_hgnc_from_up():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 phosphorylates ERK2.',
                                offline=offline)
        assert len(rp.statements) == 1
        st = rp.statements[0]
        (map2k1, mapk1) = st.agent_list()
        assert map2k1.name == 'MAP2K1'
        assert map2k1.db_refs['HGNC'] == '6840'
        assert map2k1.db_refs['UP'] == 'Q02750'
        assert mapk1.name == 'MAPK1'
        assert mapk1.db_refs['HGNC'] == '6871'
        assert mapk1.db_refs['UP'] == 'P28482'
        assert unicode_strs(rp.statements)

def assert_pmid(stmt):
    for ev in stmt.evidence:
        assert(ev.pmid is not None)
        assert(not ev.pmid.startswith('api'))
        assert(not ev.pmid.startswith('PMID'))

