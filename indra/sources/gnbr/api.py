__all__ = ['process_gene_gene']

import pandas
from .processor import GnbrGeneGeneProcessor

def process_gene_gene(part1_path, part2_path):
    df1 = pandas.read_csv(part1_path, sep='\t')
    df2 = pandas.read_csv(part2_path, sep='\t', header=None)
    gp = GnbrGeneGeneProcessor(df1, df2)
    gp.extract_activations()
    return gp