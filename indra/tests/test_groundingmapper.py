from indra.preassembler.grounding_mapper import default_mapper as gm
from indra.preassembler.grounding_mapper import GroundingMapper
from indra.preassembler.grounding_mapper.analysis import *
from indra.preassembler.grounding_mapper.gilda import ground_statements
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
    assert unicode_strs((akt, stmt, gm, mapped_akt))


def test_map_standardize_up_hgnc():
    a1 = Agent('MAPK1', db_refs={'HGNC': '6871'})
    a2 = Agent('MAPK1', db_refs={'UP': 'P28482'})
    stmt = Phosphorylation(a1, a2)
    mapped_stmts = gm.map_stmts([stmt])
    assert len(mapped_stmts) == 1
    st = mapped_stmts[0]
    assert st.enz.db_refs['HGNC'] == st.sub.db_refs['HGNC']
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
    akt_indra = Agent('pkbA', db_refs={'TEXT': 'Akt', 'FPLX': 'AKT family',
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
    assert renamed_stmts[0].sub.name == 'AKT family'
    # Then on the HGNC lookup from Uniprot
    assert renamed_stmts[1].sub.name == 'AKT1', renamed_stmts[1].sub.name
    # Don't fall back on text if there's no grounding
    assert renamed_stmts[2].sub.name == 'pkbA'
    assert renamed_stmts[3].sub.name == 'tat'
    assert unicode_strs((akt_indra, akt_hgnc_from_up, akt_other,
                         tat_up_no_hgnc, stmts, gm, renamed_stmts))


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
    assert unicode_strs((erk, stmt, gm, mapped_stmts, mapped_erk))


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
    assert unicode_strs((erk, stmt, gm, mapped_stmts, mapped_erk))


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
    assert unicode_strs((erk, stmt, gm, mapped_stmts, mapped_erk))


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
    assert ag.db_refs['HGNC'] == '377'
    assert ag.db_refs['UP'] == 'O43687'
    assert ag.name == 'AKAP7'


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
    assert unicode_strs((gag, stmt, gm, mapped_stmts, mapped_gag))


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
    assert unicode_strs((no_gn, stmt, gm, mapped_stmts, mapped_ag))


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
    assert ms.sub.db_refs == {'TEXT': 'NF-kappaB p65', 'UP': 'Q04206',
                              'HGNC': '9955'}


def test_map_agent():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    p_erk = Agent('P-ERK', db_refs={'TEXT': 'p-ERK'})
    stmt = Complex([erk, p_erk])
    mapped_stmts = gm.map_stmts([stmt])
    mapped_ag = mapped_stmts[0].members[1]
    assert mapped_ag.name == 'ERK'
    assert mapped_ag.db_refs.get('FPLX') == 'ERK'


def test_name_standardize_hgnc_up():
    a1 = Agent('x', db_refs={'HGNC': '9387'})
    GroundingMapper.standardize_agent_name(a1, True)
    assert a1.name == 'PRKAG3'
    a1 = Agent('x', db_refs={'UP': 'Q9UGI9'})
    GroundingMapper.standardize_agent_name(a1, True)
    assert a1.name == 'PRKAG3'
    a1 = Agent('x', db_refs={'UP': 'Q8BGM7'})
    GroundingMapper.standardize_agent_name(a1, True)
    assert a1.name == 'Prkag3'


def test_name_standardize_chebi():
    a1 = Agent('x', db_refs={'CHEBI': '15996'})
    GroundingMapper.standardize_agent_name(a1, False)
    assert a1.name == 'GTP'


def test_name_standardize_go():
    a1 = Agent('x', db_refs={'GO': 'GO:0006915'})
    GroundingMapper.standardize_agent_name(a1, False)
    assert a1.name == 'apoptotic process'


def test_name_standardize_mesh():
    a1 = Agent('x', db_refs={'MESH': 'D008545'})
    GroundingMapper.standardize_agent_name(a1, False)
    assert a1.name == 'Melanoma', a1.name


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
    assert mapped_stmts1[0].sub.name == 'ESR1'
    assert mapped_stmts1[0].sub.db_refs['HGNC'] == '3467'
    assert mapped_stmts1[0].sub.db_refs['UP'] == 'P03372'

    mapped_stmts2 = gm.map_stmts([stmt2])
    assert mapped_stmts2[0].obj.name == 'endoplasmic reticulum', \
        mapped_stmts2[0].obj.name
    assert mapped_stmts2[0].obj.db_refs['GO'] == 'GO:0005783', \
        mapped_stmts2[0].obj.db_refs

    annotations = mapped_stmts2[0].evidence[0].annotations
    assert 'GO:GO:0005783' in annotations['agents']['adeft'][1]


def test_misgrounding():
    baz1 = Agent('ZNF214', db_refs={'TEXT': 'baz1', 'HGNC': '13006'})
    stmt = Phosphorylation(None, baz1)
    stmts = gm.map_stmts([stmt])
    stmt = stmts[0]
    assert len(stmt.sub.db_refs) == 1, stmt.sub.db_refs
    assert stmt.sub.db_refs['TEXT'] == 'baz1'
    assert stmt.sub.name == 'baz1'


def test_ground_gilda():
    mek = Agent('Mek', db_refs={'TEXT': 'MEK'})
    erk = Agent('Erk1', db_refs={'TEXT': 'Erk1'})
    stmt = Phosphorylation(mek, erk)
    ground_statements([stmt])
    assert stmt.enz.name == 'MEK', stmt.enz
    assert stmt.enz.db_refs['FPLX'] == 'MEK'
    assert stmt.sub.name == 'MAPK3'
    assert stmt.sub.db_refs['HGNC'] == '6877'
    assert stmt.sub.db_refs['UP'] == 'P27361'
