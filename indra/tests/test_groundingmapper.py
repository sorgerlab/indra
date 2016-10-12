from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.preassembler.grounding_mapper import *
from indra.statements import Agent, Phosphorylation, Evidence
from indra.util import unicode_strs

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
    assert mapped_akt.db_refs['BE'] == 'AKT'
    assert unicode_strs((akt, stmt, gm, mapped_akt))

def test_renaming():
    akt_indra = Agent('pkbA', db_refs={'TEXT': 'Akt', 'BE':'AKT family',
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

if __name__ == '__main__':
    test_save_sentences_unicode()

