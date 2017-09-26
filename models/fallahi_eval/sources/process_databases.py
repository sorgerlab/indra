from indra.util import _require_python3
from indra.sources import biopax
from indra.tools.gene_network import GeneNetwork
import indra.tools.assemble_corpus as ac
from util import prefixed_pkl, based, basen


phosphosite_owl_file = 'sources/Kinase_substrates.owl'


def read_phosphosite_owl(fname=phosphosite_owl_file):
    bp = biopax.process_owl(fname)
    for stmt in bp.statements:
        for ev in stmt.evidence:
            ev.source_api = 'phosphosite'
            ev.epistemics = {'direct': True}
    return bp.statements


def build_prior(gene_names):
    """Build a corpus of prior Statements from PC and BEL."""
    gn = GeneNetwork(gene_names, basen)
    # Read BEL Statements
    bel_stmts = gn.get_bel_stmts(filter=False)
    ac.dump_statements(bel_stmts, prefixed_pkl('bel'))
    # Read Pathway Commons Statements
    database_filter = ['reactome', 'kegg', 'pid']
    biopax_stmts = gn.get_biopax_stmts(database_filter=database_filter)
    # Eliminate blacklisted interactions
    tmp_stmts = []
    for stmt in biopax_stmts:
        source_ids = [ev.source_id for ev in stmt.evidence]
        if set(source_ids) & set(biopax_blacklist):
            continue
        tmp_stmts.append(stmt)
    biopax_stmts = tmp_stmts
    ac.dump_statements(biopax_stmts, prefixed_pkl('biopax'))
    # Read Phosphosite Statements
    phosphosite_stmts = read_phosphosite_owl(phosphosite_owl_file)
    ac.dump_statements(phosphosite_stmts, prefixed_pkl('phosphosite'))


biopax_blacklist = \
    ['http://pathwaycommons.org/pc2/Catalysis_13953a072d6f992f2388d85f9059a475']
