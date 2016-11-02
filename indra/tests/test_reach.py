from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra import reach
from indra.reach.processor import ReachProcessor
from indra.util import unicode_strs

# Change this list to control what modes of
# reading are enabled in tests
offline_modes = [False, True]

def test_parse_site_text():
    text = ['threonine 185', 'thr 185', 'thr-185',
            'threonine residue 185', 'T185']
    assert unicode_strs(text)
    for t in text:
        residue, site = ReachProcessor._parse_site_text(t)
        assert(residue == 'T')
        assert(site == '185')
        assert unicode_strs((residue, site))

def test_parse_site_residue_only():
    text = ['serine residue', 'serine', 'a serine site']
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
