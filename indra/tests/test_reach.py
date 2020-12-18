import os
from nose.plugins.attrib import attr
from indra.sources import reach
from indra.sources.reach.processor import ReachProcessor
from indra.util import unicode_strs
from indra.statements import IncreaseAmount, DecreaseAmount, \
    Dephosphorylation, Complex, Phosphorylation, Translocation

# Change this list to control what modes of
# reading are enabled in tests
offline_modes = [True]


def test_parse_site_text():
    text = ['threonine 185', 'thr 185', 'thr-185',
            'threonine residue 185', 'T185']
    assert unicode_strs(text)
    for t in text:
        sites = ReachProcessor._parse_site_text(t)
        assert len(sites) == 1
        residue, site = sites[0]
        assert residue == 'T'
        assert site == '185'
        assert unicode_strs((residue, site))


def test_parse_site_text_number():
    t = '135'
    sites = ReachProcessor._parse_site_text(t)
    assert len(sites) == 1
    residue, site = sites[0]
    assert residue is None
    assert site == '135'
    assert unicode_strs(site)


def test_parse_site_text_number_first():
    t = '293T'
    sites = ReachProcessor._parse_site_text(t)
    assert len(sites) == 1
    residue, site = sites[0]
    assert residue == 'T'
    assert site == '293'
    assert unicode_strs((residue, site))


def test_parse_site_text_number_first_space():
    t = '293 T'
    sites = ReachProcessor._parse_site_text(t)
    assert len(sites) == 1
    residue, site = sites[0]
    assert residue == 'T'
    assert site == '293'
    assert unicode_strs((residue, site))


def test_parse_site_text_other_aa():
    t = 'A431'
    sites = ReachProcessor._parse_site_text(t)
    assert len(sites) == 1
    residue, site = sites[0]
    assert residue == 'A'
    assert site == '431'
    assert unicode_strs((residue, site))


def test_parse_site_residue_only():
    text = ['serine residue', 'serine', 'a serine site', 's', 'ser']
    assert unicode_strs(text)
    for t in text:
        sites = ReachProcessor._parse_site_text(t)
        assert len(sites) == 1
        residue, site = sites[0]
        assert unicode_strs((residue, site))
        assert residue == 'S'
        assert site is None


def test_parse_site_multiple():
    sites = ReachProcessor._parse_site_text('638/641')
    assert len(sites) == 2
    assert sites[0][0] is None
    assert sites[0][1] == '638'
    assert sites[1][0] is None
    assert sites[1][1] == '641'

    sites = ReachProcessor._parse_site_text('992,1068')
    assert len(sites) == 2
    assert sites[0][0] is None
    assert sites[0][1] == '992'
    assert sites[1][0] is None
    assert sites[1][1] == '1068'

    sites = ReachProcessor._parse_site_text('Y1221/1222')
    assert len(sites) == 2
    assert sites[0][0] == 'Y'
    assert sites[0][1] == '1221'
    assert sites[1][0] == 'Y'
    assert sites[1][1] == '1222'

    sites = ReachProcessor._parse_site_text('Tyr-577/576')
    assert len(sites) == 2
    assert sites[0][0] == 'Y'
    assert sites[0][1] == '577'
    assert sites[1][0] == 'Y'
    assert sites[1][1] == '576'

    sites = ReachProcessor._parse_site_text('S199/S202/T205')
    assert len(sites) == 3
    assert sites[0][0] == 'S'
    assert sites[0][1] == '199'
    assert sites[1][0] == 'S'
    assert sites[1][1] == '202'
    assert sites[2][0] == 'T'
    assert sites[2][1] == '205'

    sites = ReachProcessor._parse_site_text('S199/202/T205')
    assert len(sites) == 3
    assert sites[0][0] == 'S'
    assert sites[0][1] == '199'
    assert sites[1][0] is None
    assert sites[1][1] == '202'
    assert sites[2][0] == 'T'
    assert sites[2][1] == '205'

    sites = ReachProcessor._parse_site_text('S199/202/205')
    assert len(sites) == 3
    assert sites[0][0] == 'S'
    assert sites[0][1] == '199'
    assert sites[1][0] == 'S'
    assert sites[1][1] == '202'
    assert sites[2][0] == 'S'
    assert sites[2][1] == '205'


def test_phosphorylate():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 phosphorylates ERK2.', offline=offline)
        assert rp is not None
        assert len(rp.statements) == 1
        s = rp.statements[0]
        assert (s.enz.name == 'MAP2K1')
        assert (s.sub.name == 'MAPK1')
        assert unicode_strs(rp.statements)


def test_indirect_phosphorylate():
    txt = 'DUSP decreases the phosphorylation of ERK.'
    for offline in offline_modes:
        rp = reach.process_text(txt, offline=offline)
        assert rp is not None
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
        assert rp is not None
        assert len(rp.statements) == 1
        s = rp.statements[0]
        assert isinstance(s, IncreaseAmount)
        assert (s.subj.name == 'ERK')
        assert (s.obj.name == 'DUSP')
        assert unicode_strs(rp.statements)
        rp = reach.process_text('ERK decreases the amount of DUSP.',
                                offline=offline)
        assert len(rp.statements) == 1
        s = rp.statements[0]
        assert isinstance(s, DecreaseAmount)
        assert (s.subj.name == 'ERK')
        assert (s.obj.name == 'DUSP')
        assert unicode_strs(rp.statements)


def test_multiple_enzymes():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 and MEK2 phosphorylate ERK1.',
                                offline=offline)
        assert rp is not None
        assert len(rp.statements) == 2
        stmts = sorted(rp.statements, key=lambda x: x.enz.name)
        assert stmts[0].enz.name == 'MAP2K1', stmts
        assert stmts[1].enz.name == 'MAP2K2', stmts
        assert stmts[0].sub.name == 'MAPK3', stmts
        assert stmts[1].sub.name == 'MAPK3', stmts


def test_activate():
    for offline in offline_modes:
        rp = reach.process_text('HRAS activates BRAF.', offline=offline)
        assert rp is not None
        assert len(rp.statements) == 1
        s = rp.statements[0]
        assert (s.subj.name == 'HRAS')
        assert (s.obj.name == 'BRAF')
        assert unicode_strs(rp.statements)


def test_reg_amount_complex_controller():
    txt = 'The FOS-JUN complex increases the amount of ZEB2.'
    for offline in offline_modes:
        rp = reach.process_text(txt, offline=offline)
        assert rp is not None
        assert len(rp.statements) == 2
        cplx = [s for s in rp.statements if isinstance(s, Complex)][0]
        regam = [s for s in rp.statements if isinstance(s, IncreaseAmount)][0]
        assert {a.name for a in cplx.members} < {'FOS_family',
                                                 # Old version: JUN, new:
                                                 # JUN_family
                                                 'JUN_family', 'JUN'}, cplx
        assert len(regam.subj.bound_conditions) == 1
        assert unicode_strs(rp.statements)


def test_bind():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 binds ERK2.', offline=offline)
        assert rp is not None
        assert len(rp.statements) == 1
        assert unicode_strs(rp.statements)


def test_fplx_grounding():
    for offline in offline_modes:
        rp = reach.process_text('MEK activates ERK.', offline=offline)
        assert rp is not None
        assert len(rp.statements) == 1
        assert unicode_strs(rp.statements)
        if offline is True:
            st = rp.statements[0]
            assert st.subj.db_refs.get('FPLX') == 'MEK'
            assert st.obj.db_refs.get('FPLX') == 'ERK'


def test_conversions():
    here = os.path.dirname(os.path.abspath(__file__))
    test_file = os.path.join(here, 'reach_conversion.json')
    rp = reach.process_json_file(test_file)
    assert rp is not None
    assert len(rp.statements) == 1
    stmt = rp.statements[0]
    assert stmt.subj.name == 'ACE'
    assert len(stmt.obj_from) == 1
    assert stmt.obj_from[0].name == 'angiotensin-I'
    assert stmt.obj_to[0].name == 'angiotensin-II'


def test_activity():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 activates ERK2.', offline=offline)
        assert rp is not None
        assert len(rp.statements) == 1
        assert unicode_strs(rp.statements)


def test_mutation():
    for offline in offline_modes:
        rp = reach.process_text('BRAF(V600E) phosphorylates MEK.',
                                offline=offline)
        assert rp is not None
        assert len(rp.statements) == 1
        braf = rp.statements[0].enz
        assert braf.name == 'BRAF'
        assert len(braf.mutations) == 1
        assert braf.mutations[0].position == '600'
        assert braf.mutations[0].residue_from == 'V'
        assert braf.mutations[0].residue_to == 'E'
        assert unicode_strs(rp.statements)


def test_parse_mutation():
    mut = ReachProcessor._parse_mutation('V600E')
    assert mut.residue_from == 'V'
    assert mut.position == '600'
    assert mut.residue_to == 'E'

    mut = ReachProcessor._parse_mutation('Leu174Arg')
    assert mut.residue_from == 'L'
    assert mut.position == '174'
    assert mut.residue_to == 'R'

    mut = ReachProcessor._parse_mutation('val34leu')
    assert mut.residue_from == 'V'
    assert mut.position == '34'
    assert mut.residue_to == 'L'


def test_process_unicode():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 binds ERK2\U0001F4A9.', offline=offline)
        assert rp is not None
        assert unicode_strs(rp.statements)


@attr('slow')
def test_process_pmc():
    for offline in offline_modes:
        rp = reach.process_pmc('PMC4338247', offline=offline)
        assert rp is not None
        for stmt in rp.statements:
            assert_pmid(stmt)
        assert unicode_strs(rp.statements)


def test_process_unicode_abstract():
    for offline in offline_modes:
        rp = reach.process_pubmed_abstract('27749056', offline=offline)
        assert rp is not None
        assert unicode_strs(rp.statements)


def test_hgnc_from_up():
    for offline in offline_modes:
        rp = reach.process_text('MEK1 phosphorylates ERK2.',
                                offline=offline)
        assert rp is not None
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
        assert ev.pmid is not None
        assert not ev.pmid.startswith('api')
        assert not ev.pmid.startswith('PMID')


def test_process_mod_condition1():
    test_cases = [
        ('phosphorylated MEK1 activates ERK1.',
         'phosphorylation', None, None, True),
        ('MEK1 that is phosphorylated on serine activates ERK1.',
         'phosphorylation', 'S', None, True),
        ('MEK1 that is phosphorylated on S222 activates ERK1.',
         'phosphorylation', 'S', '222', True),
        ]
    for offline in offline_modes:
        for sentence, mod_type, residue, position, is_modified in test_cases:
            rp = reach.process_text(sentence, offline=offline)
            assert rp is not None
            assert len(rp.statements) == 1
            mcs = rp.statements[0].subj.mods
            assert len(mcs) == 1, 'No mods for %s' % sentence
            assert mcs[0].mod_type == mod_type, mcs
            assert mcs[0].residue == residue, mcs
            assert mcs[0].position == position, mcs
            assert mcs[0].is_modified == is_modified, mcs


def test_get_db_refs_up_human():
    entity_term = {
        'text': 'Ikaros',
        'xrefs': [{'namespace': 'uniprot', 'id': 'Q13422',
                   'object-type': 'db-reference'}]
        }
    db_refs = ReachProcessor._get_db_refs(entity_term)
    assert db_refs == {'UP': 'Q13422', 'HGNC': '13176',
                       'TEXT': 'Ikaros', 'MESH': 'C497239',
                       'EGID': '10320'}, db_refs


def test_get_db_refs_up_non_human():
    entity_term = {
        'text': 'MYC',
        'xrefs': [{'namespace': 'uniprot', 'id': 'Q9MZT7',
                   'object-type': 'db-reference'}]
        }
    db_refs = ReachProcessor._get_db_refs(entity_term)
    assert db_refs == {'UP': 'Q9MZT7', 'TEXT': 'MYC'}, db_refs


def test_get_agent_coordinates_phosphorylation():
    test_case = ('This sentence is filler. '
                 'Two filler sentences will work. '
                 'MEK that is phosphorylated phosphorylates ERK.')
    for offline in offline_modes:
        rp = reach.process_text(test_case, offline=offline)
        assert rp is not None
        stmt = rp.statements[0]
        annotations = stmt.evidence[0].annotations

        coords = [(0, 3), (42, 45)]
        assert annotations['agents']['coords'] == coords


def test_get_agent_coordinates_activation():
    test_case = 'MEK1 activates ERK2'
    for offline in offline_modes:
        rp = reach.process_text(test_case, offline=offline)
        assert rp is not None
        stmt = rp.statements[0]
        annotations = stmt.evidence[0].annotations
        coords = [(0, 4), (15, 19)]
        assert annotations['agents']['coords'] == coords


def test_get_agent_coordinates_regulate_amount():
    test_case = 'ERK increases the transcription of DUSP'
    for offline in offline_modes:
        rp = reach.process_text(test_case, offline=offline)
        assert rp is not None
        stmt = rp.statements[0]
        annotations = stmt.evidence[0].annotations
        coords = [(0, 3), (35, 39)]
        assert annotations['agents']['coords'] == coords


def test_get_agent_coordinates_binding():
    test_case = 'Everyone has observed that MEK1 binds ERK2'
    for offline in offline_modes:
        rp = reach.process_text(test_case, offline=offline)
        assert rp is not None
        stmt = rp.statements[0]
        annotations = stmt.evidence[0].annotations
        coords = [(27, 31), (38, 42)]
        assert annotations['agents']['coords'] == coords


def test_get_agent_coordinates_translocation():
    test_case = ('The length of time that ERK phosphorylation '
                 'is sustained may determine whether active ERK '
                 'translocates to the nucleus promotes cell growth.')
    for offline in offline_modes:
        rp = reach.process_text(test_case, offline=offline)
        assert rp is not None
        stmt = [stmt for stmt in rp.statements if
                isinstance(stmt, Translocation)][0]
        annotations = stmt.evidence[0].annotations
        coords = [(86, 89)]
        assert annotations['agents']['coords'] == coords


def test_get_agent_coordinates_phosphorylation_missing_controller():
    test_case = ('The ability of uPA and PAI-1 complex to induce '
                 'sustained ERK phosphorylation in MCF-7 cells '
                 'requires the recruitment of uPAR and VLDLr, which '
                 'function cooperatively')
    for offline in offline_modes:
        rp = reach.process_text(test_case, offline=offline)
        assert rp is not None
        phos_stmts = [stmt for stmt in rp.statements if
                      isinstance(stmt, Phosphorylation)]
        assert phos_stmts, rp.statements
        stmt = phos_stmts[0]
        annotations = stmt.evidence[0].annotations
        coords = [None, (57, 60)]
        assert annotations['agents']['coords'] == coords


def test_amount_embedded_in_activation():
    here = os.path.dirname(os.path.abspath(__file__))
    test_file = os.path.join(here, 'reach_act_amt.json')
    rp = reach.process_json_file(test_file)
    assert rp is not None
    assert len(rp.statements) == 1
    assert isinstance(rp.statements[0], IncreaseAmount)
    assert rp.statements[0].subj is not None
    assert rp.statements[0].obj is not None


def test_phosphorylation_regulation():
    here = os.path.dirname(os.path.abspath(__file__))
    test_file = os.path.join(here, 'reach_reg_phos.json')
    rp = reach.process_json_file(test_file)
    assert rp is not None
    assert len(rp.statements) == 1
    stmt = rp.statements[0]
    assert isinstance(stmt, Phosphorylation), stmt
    assert not stmt.sub.mods
