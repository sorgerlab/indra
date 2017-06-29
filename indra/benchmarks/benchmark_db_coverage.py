from indra import bel
from indra import biopax


def get_bel_statistics(bp):
    if bp is None:
        return 0, 0
    num_all = len(bp.all_direct_stmts) - len(bp.degenerate_stmts)
    num_extracted = len(bp.converted_direct_stmts)
    return num_all, num_extracted

def get_biopax_statistics(bp):
    if bp is None:
        return 0, 0
    num_all, num_extracted = bp.get_coverage()
    return num_all, num_extracted

signaling = ['MAPK1', 'AKT1', 'JAK1', 'MAPK8', 'CTNNB1']
genereg = ['MYC', 'TP53', 'STAT3', 'FOXO3', 'JUN']
metabolism = ['IDH1', 'PFKL', 'DHFR', 'GLUL', 'NOS1']

#all_genes = signaling + genereg + metabolism
all_genes = metabolism
stats = {'bel': {}, 'biopax': {}}
for gene in all_genes:
    print('%s\n=========' % gene)
    belp = bel.process_ndex_neighborhood([gene])
    num_all, num_extracted = get_bel_statistics(belp)
    stats['bel'][gene] = (num_all, num_extracted)
    print(num_all, num_extracted)
    biopp = biopax.process_pc_neighborhood([gene])
    num_all, num_extracted = get_biopax_statistics(biopp)
    print(num_all, num_extracted)
    stats['biopax'][gene] = (num_all, num_extracted)
