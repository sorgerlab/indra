from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import unittest
from nose.plugins.attrib import attr
from indra.sources import reach
from indra.sources.reach.processor import ReachProcessor
from indra.util import unicode_strs
from indra.statements import IncreaseAmount, DecreaseAmount, \
    Dephosphorylation, Complex

# Change this list to control what modes of
# reading are enabled in tests
offline_modes = [True]


def test_parse_site_text():
    text = ['threonine 185', 'thr 185', 'thr-185',
            'threonine residue 185', 'T185']
    assert unicode_strs(text)
    for t in text:
        sites = ReachProcessor._parse_site_text(t)
        assert(len(sites) == 1)
        residue, site = sites[0]
        assert(residue == 'T')
        assert(site == '185')
        assert unicode_strs((residue, site))


def test_parse_site_text_number():
    t = '135'
    sites = ReachProcessor._parse_site_text(t)
    assert(len(sites) == 1)
    residue, site = sites[0]
    assert(residue is None)
    assert(site == '135')
    assert(unicode_strs(site))


def test_parse_site_text_number_first():
    t = '293T'
    sites = ReachProcessor._parse_site_text(t)
    assert(len(sites) == 1)
    residue, site = sites[0]
    assert(residue == 'T')
    assert(site == '293')
    assert(unicode_strs((residue, site)))


def test_parse_site_text_number_first_space():
    t = '293 T'
    sites = ReachProcessor._parse_site_text(t)
    assert(len(sites) == 1)
    residue, site = sites[0]
    assert(residue == 'T')
    assert(site == '293')
    assert(unicode_strs((residue, site)))


def test_parse_site_text_other_aa():
    t = 'A431'
    sites = ReachProcessor._parse_site_text(t)
    assert(len(sites) == 1)
    residue, site = sites[0]
    assert(residue == 'A')
    assert(site == '431')
    assert(unicode_strs((residue, site)))


def test_parse_site_residue_only():
    text = ['serine residue', 'serine', 'a serine site', 's', 'ser']
    assert unicode_strs(text)
    for t in text:
        sites = ReachProcessor._parse_site_text(t)
        assert(len(sites) == 1)
        residue, site = sites[0]
        assert unicode_strs((residue, site))
        assert(residue == 'S')
        assert(site is None)


def test_parse_site_multiple():
    sites = ReachProcessor._parse_site_text('638/641')
    assert(len(sites) == 2)
    assert(sites[0][0] is None)
    assert(sites[0][1] == '638')
    assert(sites[1][0] is None)
    assert(sites[1][1] == '641')

    sites = ReachProcessor._parse_site_text('992,1068')
    assert(len(sites) == 2)
    assert(sites[0][0] is None)
    assert(sites[0][1] == '992')
    assert(sites[1][0] is None)
    assert(sites[1][1] == '1068')

    sites = ReachProcessor._parse_site_text('Y1221/1222')
    assert(len(sites) == 2)
    assert(sites[0][0] == 'Y')
    assert(sites[0][1] == '1221')
    assert(sites[1][0] == 'Y')
    assert(sites[1][1] == '1222')

    sites = ReachProcessor._parse_site_text('Tyr-577/576')
    assert(len(sites) == 2)
    assert(sites[0][0] == 'Y')
    assert(sites[0][1] == '577')
    assert(sites[1][0] == 'Y')
    assert(sites[1][1] == '576')

    sites = ReachProcessor._parse_site_text('S199/S202/T205')
    assert(len(sites) == 3)
    assert(sites[0][0] == 'S')
    assert(sites[0][1] == '199')
    assert(sites[1][0] == 'S')
    assert(sites[1][1] == '202')
    assert(sites[2][0] == 'T')
    assert(sites[2][1] == '205')

    sites = ReachProcessor._parse_site_text('S199/202/T205')
    assert(len(sites) == 3)
    assert(sites[0][0] == 'S')
    assert(sites[0][1] == '199')
    assert(sites[1][0] is None)
    assert(sites[1][1] == '202')
    assert(sites[2][0] == 'T')
    assert(sites[2][1] == '205')

    sites = ReachProcessor._parse_site_text('S199/202/205')
    assert(len(sites) == 3)
    assert(sites[0][0] == 'S')
    assert(sites[0][1] == '199')
    assert(sites[1][0] == 'S')
    assert(sites[1][1] == '202')
    assert(sites[2][0] == 'S')
    assert(sites[2][1] == '205')


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
        assert s.evidence[0].epistemics.get('direct') is False


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


def test_reg_amount_complex_controller():
    txt = 'The FOS-JUN complex increases the amount of ZEB2.'
    for offline in offline_modes:
        rp = reach.process_text(txt, offline=offline)
        assert len(rp.statements) == 2
        cplx = [s for s in rp.statements if isinstance(s, Complex)][0]
        regam = [s for s in rp.statements if isinstance(s, IncreaseAmount)][0]
        assert {a.name for a in cplx.members} == {'FOS_family', 'JUN'}
        assert len(regam.subj.bound_conditions) == 1
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
        if offline is True:
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


def test_parse_mutation():
    mut = ReachProcessor._parse_mutation('V600E')
    assert(mut.residue_from == 'V')
    assert(mut.position == '600')
    assert(mut.residue_to == 'E')

    mut = ReachProcessor._parse_mutation('Leu174Arg')
    assert(mut.residue_from == 'L')
    assert(mut.position == '174')
    assert(mut.residue_to == 'R')

    mut = ReachProcessor._parse_mutation('val34leu')
    assert(mut.residue_from == 'V')
    assert(mut.position == '34')
    assert(mut.residue_to == 'L')


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


def test_process_mod_condition1():
    test_cases = [
        ('MEK1 activates ERK1 that is phosphorylated.',
         'phosphorylation', None, None, True),
        ('MEK1 activates ERK1 that is phosphorylated on tyrosine.',
         'phosphorylation', 'Y', None, True),
        ('MEK1 activates ERK1 that is phosphorylated on Y185.',
         'phosphorylation', 'Y', '185', True),
        ]
    for offline in offline_modes:
        for sentence, mod_type, residue, position, is_modified in test_cases:
            rp = reach.process_text(sentence)
            assert rp is not None
            assert len(rp.statements) == 1
            mcs = rp.statements[0].obj.mods
            assert len(mcs) == 1
            assert mcs[0].mod_type == mod_type
            assert mcs[0].residue == residue
            assert mcs[0].position == position
            assert mcs[0].is_modified == is_modified


def test_get_db_refs_up_human():
    entity_term = {
        'text': 'Ikaros',
        'xrefs': [{'namespace': 'uniprot', 'id': 'Q13422',
                   'object-type': 'db-reference'}]
        }
    name, db_refs = ReachProcessor._get_db_refs(entity_term)
    assert name == 'IKZF1', name
    assert db_refs == {'UP': 'Q13422', 'HGNC': '13176',
                       'TEXT': 'Ikaros'}, db_refs


def test_get_db_refs_up_non_human():
    entity_term = {
        'text': 'MYC',
        'xrefs': [{'namespace': 'uniprot', 'id': 'Q9MZT7',
                   'object-type': 'db-reference'}]
        }
    name, db_refs = ReachProcessor._get_db_refs(entity_term)
    assert name == 'MYC', name
    assert db_refs == {'UP': 'Q9MZT7', 'TEXT': 'MYC'}, db_refs
