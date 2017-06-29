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

def print_stats(stats):
    header = ['Gene', 'BPXtot', 'BPXext', 'BPX\%', 'BELtot', 'BELext', 'BEL\%']
    print(' & '.join(header) + '\\\\')
    for gene, gene_stats in stats.items():
        bpx_vals = ['%d' % gene_stats['biopax'][0],
                    '%d' % gene_stats['biopax'][1],
                    ('%.1f\\%%' % (100.0*gene_stats['biopax'][1] /
                                   gene_stats['biopax'][0])
                                 if gene_stats['biopax'][0] else 'N/A')]
        bel_vals = ['%d' % gene_stats['bel'][0],
                    '%d' % gene_stats['bel'][1],
                    ('%.1f\\%%' % (100.0*gene_stats['bel'][1] /
                                   gene_stats['bel'][0])
                                if gene_stats['bel'][0] else 'N/A')]
        print(' & '.join([gene] + bpx_vals + bel_vals) + '\\\\')

if __name__ == '__main__':
    signaling = ['MAPK1', 'AKT1', 'JAK1', 'GNAS', 'CTNNB1']
    genereg = ['MYC', 'TP53', 'STAT3', 'FOXO3', 'JUN']
    metabolism = ['IDH1', 'PFKL', 'DHFR', 'GLUL', 'NOS1']

    all_genes = signaling + genereg + metabolism
    stats = {g: {} for g in all_genes}
    for gene in all_genes:
        print('%s\n======' % gene)
        belp = bel.process_ndex_neighborhood([gene])
        num_all, num_extracted = get_bel_statistics(belp)
        stats[gene]['bel'] = (num_all, num_extracted)
        print(num_all, num_extracted)
        biopp = biopax.process_pc_neighborhood([gene])
        num_all, num_extracted = get_biopax_statistics(biopp)
        print(num_all, num_extracted)
        stats[gene]['biopax'] = (num_all, num_extracted)

    print_stats(stats)
