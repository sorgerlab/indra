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
    gm.map_agents([stmt])



