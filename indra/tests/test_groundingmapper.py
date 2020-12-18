from indra.preassembler.grounding_mapper import default_mapper as gm
from indra.preassembler.grounding_mapper import GroundingMapper
from indra.preassembler.grounding_mapper.analysis import *
from indra.preassembler.grounding_mapper.gilda import ground_statements, \
    get_gilda_models, ground_statement
from indra.statements import Agent, Phosphorylation, Complex, Inhibition, \
    Evidence, BoundCondition
from indra.util import unicode_strs
from nose.tools import raises
from nose.plugins.attrib import attr


def test_simple_mapping():
    akt = Agent('pkbA', db_refs={'TEXT': 'Akt', 'UP': 'XXXXXX'})
    stmt = Phosphorylation(None, akt)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    mapped_akt = mapped_stmts[0].sub
    assert mapped_akt.db_refs['TEXT'] == 'Akt'
    assert mapped_akt.db_refs['FPLX'] == 'AKT'


def test_map_standardize_up_hgnc():
    a1 = Agent('MAPK1', db_refs={'HGNC': '6871'})
    a2 = Agent('MAPK1', db_refs={'UP': 'P28482'})
    stmt = Phosphorylation(a1, a2)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    st = mapped_stmts[0]
    assert st.enz.db_refs['HGNC'] == st.sub.db_refs['HGNC'], \
        (st.enz.db_refs, st.sub.db_refs)
    assert st.enz.db_refs['UP'] == st.sub.db_refs['UP']


def test_map_standardize_chebi_pc():
    a1 = Agent('X', db_refs={'PUBCHEM': '42611257'})
    a2 = Agent('Y', db_refs={'CHEBI': 'CHEBI:63637'})
    stmt = Phosphorylation(a1, a2)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    st = mapped_stmts[0]
    assert st.enz.db_refs['PUBCHEM'] == st.sub.db_refs['PUBCHEM'], \
        (st.enz.db_refs, st.sub.db_refs)
    assert st.enz.db_refs['CHEBI'] == st.sub.db_refs['CHEBI'], \
        (st.enz.db_refs, st.sub.db_refs)
    assert st.enz.name == 'vemurafenib'
    assert st.sub.name == 'vemurafenib'


def test_map_standardize_chebi_hmdb():
    a1 = Agent('X', db_refs={'HMDB': 'HMDB0000122'})
    a2 = Agent('Y', db_refs={'CHEBI': 'CHEBI:4167'})
    stmt = Phosphorylation(a1, a2)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    st = mapped_stmts[0]
    assert st.enz.db_refs['CHEBI'] == st.sub.db_refs['CHEBI'], \
        (st.enz.db_refs, st.sub.db_refs)
    assert st.enz.name == 'D-glucopyranose', st.enz
    assert st.sub.name == 'D-glucopyranose', st.sub


def test_bound_condition_mapping():
    # Verify that the grounding mapper grounds the agents within a bound
    # condition
    akt = Agent('pkbA', db_refs={'TEXT': 'Akt', 'UP': 'XXXXXX'})
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    akt.bound_conditions = [BoundCondition(erk)]

    stmt = Phosphorylation(None, akt)

    mapped_stmts = gm.map_stmts([stmt])

    s = mapped_stmts[0]
    mapped_akt = mapped_stmts[0].sub
    mapped_erk = mapped_akt.bound_conditions[0].agent

    assert mapped_akt.db_refs['TEXT'] == 'Akt'
    assert mapped_akt.db_refs['FPLX'] == 'AKT'

    assert mapped_erk.db_refs['TEXT'] == 'ERK1'
    assert mapped_erk.db_refs['HGNC'] == '6877'
    assert mapped_erk.db_refs['UP'] == 'P27361'


def test_bound_condition_mapping_multi():
    # Test with multiple agents
    akt = Agent('pkbA', db_refs={'TEXT': 'Akt', 'UP': 'XXXXXX'})
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    akt.bound_conditions = [BoundCondition(erk)]
    stmt = Phosphorylation(akt, erk)
    mapped_stmts = gm.map_stmts([stmt])
    s = mapped_stmts[0]
    mapped_akt = mapped_stmts[0].enz
    mapped_erk1 = mapped_akt.bound_conditions[0].agent
    mapped_erk2 = mapped_stmts[0].sub

    assert mapped_akt.db_refs['TEXT'] == 'Akt'
    assert mapped_akt.db_refs['FPLX'] == 'AKT'

    for e in (mapped_erk1, mapped_erk2):
        assert e.db_refs['TEXT'] == 'ERK1'
        assert e.db_refs['HGNC'] == '6877'
        assert e.db_refs['UP'] == 'P27361'


def test_bound_condition_mapping_agent_json():
    # Test with agent/json mapping
    akt = Agent('pkbA', db_refs={'TEXT': 'p-Akt', 'UP': 'XXXXXX'})
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    akt.bound_conditions = [BoundCondition(erk)]
    stmt = Phosphorylation(None, akt)

    mapped_stmts = gm.map_stmts([stmt])

    s = mapped_stmts[0]
    mapped_akt = mapped_stmts[0].sub
    mapped_erk = mapped_akt.bound_conditions[0].agent

    #assert mapped_akt.db_refs['TEXT'] == 'p-AKT', mapped_akt.db_refs
    assert mapped_akt.db_refs['FPLX'] == 'AKT', mapped_akt.db_refs

    assert mapped_erk.db_refs['TEXT'] == 'ERK1'
    assert mapped_erk.db_refs['HGNC'] == '6877'
    assert mapped_erk.db_refs['UP'] == 'P27361'


def test_ignore():
    agent = Agent('FA', db_refs={'TEXT': 'FA'})
    stmt = Phosphorylation(None, agent)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 0


def test_renaming():
    akt_indra = Agent('pkbA', db_refs={'TEXT': 'Akt', 'FPLX': 'AKT',
                                       'UP': 'P31749'})
    akt_hgnc_from_up = Agent('pkbA', db_refs={'TEXT': 'Akt', 'UP': 'P31749'})
    akt_other = Agent('pkbA', db_refs={'TEXT': 'Akt'})
    tat_up_no_hgnc = Agent('foo', db_refs={'TEXT': 'bar', 'UP': 'P04608'})
    stmts = [Phosphorylation(None, akt_indra),
             Phosphorylation(None, akt_hgnc_from_up),
             Phosphorylation(None, akt_other),
             Phosphorylation(None, tat_up_no_hgnc), ]
    renamed_stmts = gm.rename_agents(stmts)
    assert len(renamed_stmts) == 4
    # Should draw on BE first
    assert renamed_stmts[0].sub.name == 'AKT'
    # Then on the HGNC lookup from Uniprot
    assert renamed_stmts[1].sub.name == 'AKT1', renamed_stmts[1].sub.name
    # Don't fall back on text if there's no grounding
    assert renamed_stmts[2].sub.name == 'pkbA'
    assert renamed_stmts[3].sub.name == 'tat'


def test_save_sentences_unicode():
    mek = Agent('MEK', db_refs={'TEXT': 'MAP2K1'})
    ev = Evidence(source_api='reach', pmid='PMID000asdf',
                  text='foo\U0001F4A9bar')
    st = Phosphorylation(None, mek, evidence=[ev])
    sent = get_sentences_for_agent('MAP2K1', [st])
    assert unicode_strs(sent)
    twg = agent_texts_with_grounding([st])
    save_sentences(twg, [st], 'test_save_sentences.csv')


def test_hgnc_sym_but_not_up():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'HGNC': '6871'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    mapped_erk = mapped_stmts[0].sub
    assert mapped_erk.name == 'MAPK1'
    assert mapped_erk.db_refs['TEXT'] == 'ERK1'
    assert mapped_erk.db_refs['HGNC'] == '6871'
    assert mapped_erk.db_refs['UP'] == 'P28482'


def test_up_but_not_hgnc():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'P28482'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    mapped_erk = mapped_stmts[0].sub
    assert mapped_erk.name == 'MAPK1'
    assert mapped_erk.db_refs['TEXT'] == 'ERK1'
    assert mapped_erk.db_refs['HGNC'] == '6871'
    assert mapped_erk.db_refs['UP'] == 'P28482'


def test_hgnc_but_not_up():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'HGNC': '6871'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    mapped_erk = mapped_stmts[0].sub
    assert mapped_erk.name == 'MAPK1'
    assert mapped_erk.db_refs['TEXT'] == 'ERK1'
    assert mapped_erk.db_refs['HGNC'] == '6871'
    assert mapped_erk.db_refs['UP'] == 'P28482'


@raises(ValueError)
def test_hgnc_sym_with_no_id():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'HGNC': 'foobar'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])


@raises(ValueError)
def test_up_and_invalid_hgnc_sym():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'P28482', 'HGNC': 'foobar'}}
    gm = GroundingMapper(g_map)


def test_up_with_no_gene_name_with_hgnc_sym():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'A0K5Q6', 'HGNC': '6871'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    assert mapped_stmts[0].sub.db_refs['HGNC'] == '6871', \
        mapped_stmts[0].sub.db_refs
    assert mapped_stmts[0].sub.db_refs['UP'] == 'P28482', \
        mapped_stmts[0].sub.db_refs


def test_multiple_mapped_up():
    ag = Agent('xx', db_refs={'HGNC': '377', 'UP': 'O43687'})
    gm.standardize_agent_name(ag, True)
    assert ag.db_refs['HGNC'] == '377', ag.db_refs
    assert ag.db_refs['UP'] == 'O43687', ag.db_refs
    assert ag.name == 'AKAP7', ag.name


def test_up_and_mismatched_hgnc():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'P28482', 'HGNC': '6877'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    assert mapped_stmts[0].sub.db_refs['HGNC'] == '6877', \
        mapped_stmts[0].sub.db_refs
    assert mapped_stmts[0].sub.db_refs['UP'] == 'P27361', \
        mapped_stmts[0].sub.db_refs


def test_up_id_with_no_hgnc_id():
    """Non human protein"""
    gag = Agent('Gag', db_refs={'TEXT': 'Gag'})
    stmt = Phosphorylation(None, gag)
    g_map = {'Gag': {'TEXT': 'Gag', 'UP': 'P04585'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    mapped_gag = mapped_stmts[0].sub
    assert mapped_gag.name == 'gag-pol'
    assert mapped_gag.db_refs['TEXT'] == 'Gag'
    assert mapped_gag.db_refs.get('HGNC') is None
    assert mapped_gag.db_refs['UP'] == 'P04585'


def test_up_id_with_no_gene_name():
    """Expect no HGNC entry; no error raised."""
    no_gn = Agent('NoGNname', db_refs={'TEXT': 'NoGN'})
    stmt = Phosphorylation(None, no_gn)
    g_map = {'NoGN': {'TEXT': 'NoGN', 'UP': 'A0K5Q6'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    mapped_ag = mapped_stmts[0].sub
    assert mapped_ag.name == 'NoGNname'
    assert mapped_ag.db_refs['TEXT'] == 'NoGN'
    assert mapped_ag.db_refs.get('HGNC') is None
    assert mapped_ag.db_refs['UP'] == 'A0K5Q6'


def test_in_place_overwrite_of_gm():
    """Make sure HGNC lookups don't modify the original grounding map by adding
    keys."""
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'P28482'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    gmap_after_mapping = gm.grounding_map
    assert set(gmap_after_mapping['ERK1'].keys()) == set(['TEXT', 'UP'])


def test_map_entry_hgnc_and_up():
    """Make sure that HGNC symbol is replaced with HGNC ID when grounding map
    includes both UP ID and HGNC symbol."""
    rela = Agent('NF-kappaB p65', db_refs={'TEXT': 'NF-kappaB p65'})
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(erk, rela)
    g_map = {'NF-kappaB p65': {'TEXT': 'NF-kappaB p65', 'UP': 'Q04206',
                               'HGNC': '9955'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    ms = mapped_stmts[0]
    assert ms.sub.db_refs == \
           {'TEXT': 'NF-kappaB p65', 'UP': 'Q04206',
            'HGNC': '9955', 'MESH': 'D051996',
            'EGID': '5970'}, ms.sub.db_refs


def test_map_agent():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    p_erk = Agent('P-ERK', db_refs={'TEXT': 'p-ERK'})
    stmt = Complex([erk, p_erk])
    mapped_stmts = gm.map_stmts([stmt])
    mapped_ag = mapped_stmts[0].members[1]
    assert mapped_ag.name == 'ERK'
    assert mapped_ag.db_refs.get('FPLX') == 'ERK'


@attr('nonpublic')
def test_adeft_mapping():
    er1 = Agent('ER', db_refs={'TEXT': 'ER'})
    pmid1 = '30775882'
    stmt1 = Phosphorylation(None, er1, evidence=[Evidence(pmid=pmid1,
                                                          text_refs={'PMID':
                                                                     pmid1})])

    er2 = Agent('ER', db_refs={'TEXT': 'ER'})
    pmid2 = '28369137'
    stmt2 = Inhibition(None, er2, evidence=[Evidence(pmid=pmid2,
                                                     text_refs={'PMID':
                                                                pmid2})])

    mapped_stmts1 = gm.map_stmts([stmt1])
    assert mapped_stmts1[0].sub.name == 'ESR', \
        mapped_stmts1[0].sub.name
    assert mapped_stmts1[0].sub.db_refs['FPLX'] == 'ESR', \
        mapped_stmts1[0].sub.db_refs

    mapped_stmts2 = gm.map_stmts([stmt2])
    assert mapped_stmts2[0].obj.name == 'endoplasmic reticulum', \
        mapped_stmts2[0].obj.name
    assert mapped_stmts2[0].obj.db_refs['GO'] == 'GO:0005783', \
        mapped_stmts2[0].obj.db_refs

    annotations = mapped_stmts2[0].evidence[0].annotations
    assert 'GO:GO:0005783' in annotations['agents']['adeft'][1]


def test_adeft_mapping_non_pos():
    er = Agent('ER', db_refs={'TEXT': 'ER'})
    # This is an exact definition of a pos_label entry so we
    # expect that it will be applied as a grounding even though the
    # Adeft model has low precision for this label.
    ev = Evidence(text='estradiol (ER)')
    stmt = Phosphorylation(None, er, evidence=[ev])
    mapped_stmt = gm.map_stmts([stmt])[0]
    assert 'CHEBI' in mapped_stmt.sub.db_refs, mapped_stmt.evidence
    # This one is not an exact definition so we expect the grounding to
    # be stripped out.
    ev = Evidence(text='Estradiol is one of the three estrogen hormones'
                  'naturally produced in the body.')
    stmt = Phosphorylation(None, er, evidence=[ev])
    mapped_stmt = gm.map_stmts([stmt])[0]
    assert 'CHEBI' not in mapped_stmt.sub.db_refs, mapped_stmt.evidence
    # This is a non-positive label, and we expect it to be stripped out
    # whether it's an exact definition or not.
    pcs = Agent('PCS', db_refs={'TEXT': 'PCS', 'MESH': 'xxx'})
    ev = Evidence(text='physical component summary (PCS)')
    stmt = Phosphorylation(None, pcs, evidence=[ev])
    mapped_stmt = gm.map_stmts([stmt])[0]
    assert 'MESH' not in mapped_stmt.sub.db_refs, \
        (mapped_stmt.sub.db_refs, mapped_stmt.evidence)
    ev = Evidence(text='physical component summary')
    stmt = Phosphorylation(None, pcs, evidence=[ev])
    mapped_stmt = gm.map_stmts([stmt])[0]
    assert 'MESH' not in mapped_stmt.sub.db_refs, \
        (mapped_stmt.sub.db_refs, mapped_stmt.evidence)


def test_misgrounding():
    baz1 = Agent('ZNF214', db_refs={'TEXT': 'baz1', 'HGNC': '13006'})
    stmt = Phosphorylation(None, baz1)
    stmts = gm.map_stmts([stmt])
    stmt = stmts[0]
    assert len(stmt.sub.db_refs) == 1, stmt.sub.db_refs
    assert stmt.sub.db_refs['TEXT'] == 'baz1'
    assert stmt.sub.name == 'baz1'


def test_ground_gilda():
    for mode in ['web', 'local']:
        mek = Agent('Mek', db_refs={'TEXT': 'MEK'})
        erk = Agent('Erk1', db_refs={'TEXT': 'Erk1'})
        stmt = Phosphorylation(mek, erk)
        grounded_stmts = ground_statements([stmt], mode=mode)
        stmt = grounded_stmts[0]
        assert stmt.enz.name == 'MEK', stmt.enz
        assert stmt.enz.db_refs['FPLX'] == 'MEK'
        assert stmt.sub.name == 'MAPK3'
        assert stmt.sub.db_refs['HGNC'] == '6877'
        assert stmt.sub.db_refs['UP'] == 'P27361'


def test_ground_gilda_source():
    ev1 = Evidence(source_api='reach')
    ev2 = Evidence(source_api='sparser')
    ev3 = Evidence(source_api='trips')
    stmts = [Phosphorylation(None, Agent('x', db_refs={'TEXT': 'kras'}),
                             evidence=ev)
             for ev in (ev1, ev2, ev3)]
    grounded_stmts = ground_statements(stmts, sources=['trips'])
    assert grounded_stmts[0].sub.name == 'x', stmts[0]
    assert grounded_stmts[1].sub.name == 'x'
    assert grounded_stmts[2].sub.name == 'KRAS'
    grounded_stmts = ground_statements(stmts, sources=['reach', 'sparser'])
    assert all(stmt.sub.name == 'KRAS'
               for stmt in grounded_stmts[:2])


def test_gilda_ground_ungrounded():
    ag1 = Agent('x', db_refs={'TEXT': 'RAS', 'FPLX': 'RAS'})
    ag2 = Agent('x', db_refs={'TEXT': 'RAS'})
    ag3 = Agent('x', db_refs={'TEXT': 'RAS', 'XXXXX': 'XXXX'})
    stmts = [Phosphorylation(None, ag) for ag in (ag1, ag2, ag3)]
    ground_statement(stmts[0], ungrounded_only=True)
    assert ag1.name == 'x'
    ground_statement(stmts[0], ungrounded_only=False)
    assert ag1.name == 'RAS', ag1
    ground_statement(stmts[1], ungrounded_only=True)
    assert ag2.name == 'RAS'
    grounded_stmts = ground_statements([stmts[2]], ungrounded_only=True)
    assert grounded_stmts[0].sub.name == 'RAS'


def test_get_gilda_models():
    models = get_gilda_models()
    assert 'NDR1' in models


@attr('nonpublic')
def test_gilda_disambiguation():
    gm.gilda_mode = 'web'
    er1 = Agent('NDR1', db_refs={'TEXT': 'NDR1'})
    pmid1 = '18362890'
    stmt1 = Phosphorylation(None, er1,
                            evidence=[Evidence(pmid=pmid1,
                                               text_refs={'PMID': pmid1})])

    er2 = Agent('NDR1', db_refs={'TEXT': 'NDR1'})
    pmid2 = '16832411'
    stmt2 = Inhibition(None, er2,
                       evidence=[Evidence(pmid=pmid2,
                                          text_refs={'PMID': pmid2})])

    mapped_stmts1 = gm.map_stmts([stmt1])
    assert mapped_stmts1[0].sub.name == 'STK38', mapped_stmts1[0].sub.name
    assert mapped_stmts1[0].sub.db_refs['HGNC'] == '17847', \
        mapped_stmts1[0].sub.db_refs
    assert mapped_stmts1[0].sub.db_refs['UP'] == 'Q15208', \
        mapped_stmts1[0].sub.db_refs

    mapped_stmts2 = gm.map_stmts([stmt2])
    assert mapped_stmts2[0].obj.name == 'NDRG1', \
        mapped_stmts2[0].obj.name
    assert mapped_stmts2[0].obj.db_refs['HGNC'] == '7679', \
        mapped_stmts2[0].obj.db_refs
    assert mapped_stmts2[0].obj.db_refs['UP'] == 'Q92597', \
        mapped_stmts2[0].obj.db_refs

    annotations = mapped_stmts2[0].evidence[0].annotations
    assert len(annotations['agents']['gilda'][1]) == 2, \
        annotations
    assert annotations['agents']['gilda'][0] is None
    assert annotations['agents']['gilda'][1] is not None


@attr('nonpublic')
def test_gilda_disambiguation_local():
    gm.gilda_mode = 'local'
    er1 = Agent('NDR1', db_refs={'TEXT': 'NDR1'})
    pmid1 = '18362890'
    stmt1 = Phosphorylation(None, er1,
                            evidence=[Evidence(pmid=pmid1,
                                               text_refs={'PMID': pmid1})])
    mapped_stmts1 = gm.map_stmts([stmt1])
    annotations = mapped_stmts1[0].evidence[0].annotations
    assert annotations['agents']['gilda'][0] is None
    assert annotations['agents']['gilda'][1] is not None
    assert len(annotations['agents']['gilda'][1]) == 2, \
        annotations
    # This is to make sure the to_json of the ScoredMatches works
    assert annotations['agents']['gilda'][1][0]['term']['db'] == 'HGNC'


def test_grounding_map_gilda_priority():
    gm.gilda_mode = 'web'
    fetal_bovine_serum = Agent('FBS', db_refs={'TEXT': 'FBS'})
    pmid = '28536624'
    stmt = Phosphorylation(None, fetal_bovine_serum,
                           evidence=[Evidence(pmid=pmid,
                                              text_refs={'PMID': pmid})])
    mapped_stmts = gm.map_stmts([stmt])
    annotations = mapped_stmts[0].evidence[0].annotations
    # agents should not be in annotations if gilda is run. Second condition
    # added as future proofing in case some future change causes this mapping
    # to add agent annotations in the future.
    assert 'agents' not in annotations or \
        'gilda' not in annotations['agents']


def test_text_and_norm_text():
    gm.gilda_mode = 'local'

    # We should filter out ignores in both TEXT and TEXT_NORM
    ag = Agent('x', db_refs={'TEXT': 'XREF_BIBR', 'TEXT_NORM': 'ERK'})
    stmt = Phosphorylation(None, ag)
    res = gm.map_stmts([stmt])
    assert not res
    ag = Agent('x', db_refs={'TEXT': 'ERK', 'TEXT_NORM': 'XREF_BIBR'})
    stmt = Phosphorylation(None, ag)
    res = gm.map_stmts([stmt])
    assert not res

    # We should disambiguate based on both TEXT and TEXT_NORM
    ag = Agent('x', db_refs={'TEXT': 'AA', 'TEXT_NORM': 'XXX'},)
    stmt = Phosphorylation(None, ag,
                           evidence=Evidence(text='Arachidonic acid (AA)'))
    res = gm.map_stmts([stmt])
    assert res[0].sub.name == 'arachidonic acid', res[0]
    ag = Agent('x', db_refs={'TEXT': 'XXX', 'TEXT_NORM': 'AA'})
    stmt = Phosphorylation(None, ag,
                           evidence=Evidence(text='Arachidonic acid (AA)'))
    res = gm.map_stmts([stmt])
    assert res[0].sub.name == 'arachidonic acid', res[0]

    ag = Agent('x', db_refs={'TEXT': 'XXX', 'TEXT_NORM': 'ERK'})
    stmt = Phosphorylation(None, ag)
    res = gm.map_stmts([stmt])
    assert res[0].sub.name == 'ERK', res[0]

    ag = Agent('x', db_refs={'TEXT': 'ERK', 'TEXT_NORM': 'XXX'})
    stmt = Phosphorylation(None, ag)
    res = gm.map_stmts([stmt])
    assert res[0].sub.name == 'ERK', res[0]


def test_none_text_corner_case():
    ag = Agent('x', db_refs={'TEXT': None, 'TEXT_NORM': None})
    stmt = Phosphorylation(None, ag)
    res = gm.map_stmts([stmt])
    assert res[0].sub.name == 'x', res[0]
