import pickle
from indra.util import _require_python3
from indra.sources import biopax
from indra.tools.gene_network import GeneNetwork


phosphosite_owl_file = 'Kinase_substrates.owl'


def read_phosphosite_owl(fname=phosphosite_owl_file):
    bp = biopax.process_owl(fname)
    for stmt in bp.statements:
        for ev in stmt.evidence:
            ev.source_api = 'phosphosite'
            ev.epistemics = {'direct': True}
    return bp.statements


def grouped_biopax_query(gene_names, database_filter, block_size=60):
    gene_blocks = [gene_names[i:i+block_size] for i in
                   range(0, len(gene_names), block_size)]
    stmts = []
    # Run pathsfromto between pairs of blocks and pathsbetween
    # within each block. This breaks up a single call with N genes into
    # (N/block_size)*(N/blocksize) calls with block_size genes
    for genes1, genes2 in itertools.product(gene_blocks, repeat=2):
        if genes1 == genes2:
            bp = biopax.process_pc_pathsbetween(genes1,
                                            database_filter=database_filter)
        else:
            bp = biopax.process_pc_pathsfromto(genes1, genes2,
                                           database_filter=database_filter)
        stmts += bp.statements
    # Filter out blacklist
    final_stmts = []
    for stmt in stmts:
        source_ids = [ev.source_id for ev in stmt.evidence]
        if set(ev.source_id) & set(biopax_blacklist):
            continue
        final_stmts.append(stmt)
    return final_stmts


def build_prior(gene_names):
    """Build a corpus of prior Statements from PC and BEL."""
    gn = GeneNetwork(gene_names, basen)
    bel_stmts = gn.get_bel_stmts(filter=False)
    ac.dump_statements(bel_stmts, prefixed_pkl('bel'))
    # This call results in 504 error currently
    #biopax_stmts = gn.get_biopax_stmts(filter=False)
    database_filter = ['reactome', , 'kegg', 'pid']
    biopax_stmts = grouped_biopax_query(gene_names, database_filter)
    ac.dump_statements(biopax_stmts, prefixed_pkl('biopax'))
    phosphosite_stmts = read_phosphosite_owl(phosphosite_owl_file)
    ac.dump_statements(phosphosite_stmts, prefixed_pkl('phosphosite'))


biopax_blacklist = \
    ['http://pathwaycommons.org/pc2/Catalysis_13953a072d6f992f2388d85f9059a475']
