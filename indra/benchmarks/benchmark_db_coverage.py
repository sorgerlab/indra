import json
from indra.sources import bel, biopax
from collections import OrderedDict


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
    header = ['Gene', 'BPXtot', 'BPXext', 'BPX\\%', 'BELtot', 'BELext', 'BEL\\%']
    print(' & '.join(header) + '\\\\')
    for group, gene_stats in stats.items():
        for gene, gene_stats in gene_stats.items():
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
        print('\\hline')


if __name__ == '__main__':
    all_genes = \
        OrderedDict({'signaling': ['MAPK1', 'AKT1', 'JAK1', 'GNAS', 'CTNNB1'],
                     'genereg': ['MYC', 'TP53', 'STAT3', 'FOXO3', 'JUN'],
                     'metabolism': ['IDH1', 'PFKL', 'DHFR', 'GLUL', 'NOS1',
                                    'CHEBI:20506', 'CHEBI:28300', 'CHEBI:16084',
                                    'CHEBI:32816', 'CHEBI:16480']})

    stats = {group: {g: {} for g in genes}
             for group, genes in all_genes.items()}
    for group, genes in all_genes.items():
        for gene in genes:
            print('%s\n======' % gene)
            belp = bel.process_ndex_neighborhood([gene])
            num_all, num_extracted = get_bel_statistics(belp)
            stats[group][gene]['bel'] = (num_all, num_extracted)
            print(num_all, num_extracted)
            if gene.startswith('CHEBI:'):
                bpx_query = 'http://identifiers.org/chebi/' + gene
            else:
                bpx_query = gene
            biopp = biopax.process_pc_neighborhood([bpx_query])
            num_all, num_extracted = get_biopax_statistics(biopp)
            print(num_all, num_extracted)
            stats[group][gene]['biopax'] = (num_all, num_extracted)

    with open('db_coverage_stats.json', 'w') as fh:
        json.dump(stats, fh)

    print_stats(stats)
