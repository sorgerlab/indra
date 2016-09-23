from __future__ import print_function, unicode_literals
from indra.preassembler.grounding_mapper import GroundingMapper, \
                                                default_grounding_map
from indra.statements import Agent, Phosphorylation

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

