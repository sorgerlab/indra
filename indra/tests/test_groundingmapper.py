from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.preassembler.grounding_mapper import *
from indra.statements import Agent, Phosphorylation, Complex, Evidence
from indra.util import unicode_strs
from nose.tools import raises

# The grounding map
# Format:
#   - text string
#   - ID
#   - species (optional)
#   - namespace: uniprot, mesh, go, pfam, etc. (MIRIAM namespace names)
#   - type: BioProcess, Gene_or_gene_product, Family, Simple_chemical

gm = [['Akt1', 'P31749', None, 'uniprot', 'Gene_or_gene_product']]


def test_simple_mapping():
    akt = Agent('pkbA', db_refs={'TEXT': 'Akt', 'UP':'XXXXXX'})
    stmt = Phosphorylation(None, akt)
    gm = GroundingMapper(default_grounding_map)
    mapped_stmts = gm.map_agents([stmt])
    assert len(mapped_stmts) == 1
    mapped_akt = mapped_stmts[0].sub
    assert mapped_akt.db_refs['TEXT'] == 'Akt'
    assert mapped_akt.db_refs['FPLX'] == 'AKT'
    assert unicode_strs((akt, stmt, gm, mapped_akt))

def test_ignore():
    agent = Agent('FA', db_refs={'TEXT': 'FA'})
    stmt = Phosphorylation(None, agent)
    gm = GroundingMapper(default_grounding_map)
    mapped_stmts = gm.map_agents([stmt])
    assert len(mapped_stmts) == 0

def test_renaming():
    akt_indra = Agent('pkbA', db_refs={'TEXT': 'Akt', 'FPLX':'AKT family',
                                        'UP': 'P31749'})
    akt_hgnc_from_up = Agent('pkbA', db_refs={'TEXT': 'Akt', 'UP':'P31749'})
    akt_other = Agent('pkbA', db_refs={'TEXT': 'Akt'})
    tat_up_no_hgnc = Agent('foo', db_refs={'TEXT': 'bar', 'UP':'P04608'})
    stmts = [Phosphorylation(None, akt_indra),
             Phosphorylation(None, akt_hgnc_from_up),
             Phosphorylation(None, akt_other),
             Phosphorylation(None, tat_up_no_hgnc), ]
    gm = GroundingMapper(default_grounding_map)
    renamed_stmts = gm.rename_agents(stmts)
    assert len(renamed_stmts) == 4
    # Should draw on BE first
    assert renamed_stmts[0].sub.name == 'AKT family'
    # Then on the HGNC lookup from Uniprot
    assert renamed_stmts[1].sub.name == 'AKT1'
    # Don't fall back on text if there's no grounding
    assert renamed_stmts[2].sub.name == 'pkbA'
    assert renamed_stmts[3].sub.name == 'tat'
    assert unicode_strs((akt_indra, akt_hgnc_from_up, akt_other,
                         tat_up_no_hgnc, stmts, gm, renamed_stmts))

def test_save_sentences_unicode():
    mek = Agent('MEK', db_refs={'TEXT':'MAP2K1'})
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
    g_map = {'ERK1': {'TEXT': 'ERK1', 'HGNC': 'MAPK1'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])
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
    mapped_stmts = gm.map_agents([stmt])
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
    g_map = {'ERK1': {'TEXT': 'ERK1', 'HGNC': 'MAPK1'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])
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
    mapped_stmts = gm.map_agents([stmt])

@raises(ValueError)
def test_up_and_invalid_hgnc_sym():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'P28482', 'HGNC':'foobar'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])

@raises(ValueError)
def test_up_with_no_gene_name_with_hgnc_sym():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'A0K5Q6', 'HGNC':'MAPK1'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])

@raises(ValueError)
def test_up_and_mismatched_hgnc():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'P28482', 'HGNC':'MAPK3'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])

def up_id_with_no_hgnc_id():
    """Non human protein"""
    gag = Agent('Gag', db_refs={'TEXT': 'Gag'})
    stmt = Phosphorylation(None, gag)
    g_map = {'Gag': {'TEXT': 'Gag', 'UP': 'P04585'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])
    assert len(mapped_stmts) == 1
    mapped_gag = mapped_stmts[0].sub
    assert mapped_gag.name == 'gag-pol'
    assert mapped_gag.db_refs['TEXT'] == 'Gag'
    assert mapped_gag.db_refs.get('HGNC') == None
    assert mapped_gag.db_refs['UP'] == 'P04585'
    assert unicode_strs((gag, stmt, gm, mapped_stmts, mapped_gag))

def up_id_with_no_gene_name():
    """Expect no HGNC entry; no error raised."""
    no_gn = Agent('NoGNname', db_refs={'TEXT': 'NoGN'})
    stmt = Phosphorylation(None, no_gn)
    g_map = {'NoGN': {'TEXT': 'NoGN', 'UP': 'A0K5Q6'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])
    assert len(mapped_stmts) == 1
    mapped_ag = mapped_stmts[0].sub
    assert mapped_ag.name == 'NoGNname'
    assert mapped_ag.db_refs['TEXT'] == 'NoGN'
    assert mapped_ag.db_refs.get('HGNC') == None
    assert mapped_ag.db_refs['UP'] == 'A0K5Q6'
    assert unicode_strs((erk, stmt, gm, mapped_stmts, mapped_erk))

def test_in_place_overwrite_of_gm():
    """Make sure HGNC lookups don't modify the original grounding map by adding
    keys."""
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(None, erk)
    g_map = {'ERK1': {'TEXT': 'ERK1', 'UP': 'P28482'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])
    gmap_after_mapping = gm.gm
    assert set(gmap_after_mapping['ERK1'].keys()) == set(['TEXT', 'UP'])

def test_map_entry_hgnc_and_up():
    """Make sure that HGNC symbol is replaced with HGNC ID when grounding map
    includes both UP ID and HGNC symbol."""
    rela = Agent('NF-kappaB p65', db_refs={'TEXT': 'NF-kappaB p65'})
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    stmt = Phosphorylation(erk, rela)
    g_map = {'NF-kappaB p65': {'TEXT': 'NF-kappaB p65', 'UP':'Q04206',
                               'HGNC': 'RELA'}}
    gm = GroundingMapper(g_map)
    mapped_stmts = gm.map_agents([stmt])
    assert len(mapped_stmts) == 1
    ms = mapped_stmts[0]
    assert ms.sub.db_refs == {'TEXT': 'NF-kappaB p65', 'UP': 'Q04206',
                             'HGNC': '9955'}

def test_map_agent():
    erk = Agent('ERK1', db_refs={'TEXT': 'ERK1'})
    p_erk = Agent('P-ERK', db_refs={'TEXT': 'p-ERK'})
    stmt = Complex([erk, p_erk])
    gm = GroundingMapper(default_grounding_map, default_agent_map)
    mapped_stmts = gm.map_agents([stmt])
    mapped_ag = mapped_stmts[0].members[1]
    assert mapped_ag.name == 'ERK'
    assert mapped_ag.db_refs.get('FPLX') == 'ERK'

if __name__ == '__main__':
    test_map_entry_hgnc_and_up()

