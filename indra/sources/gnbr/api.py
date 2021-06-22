__all__ = ['process_gene_gene', 'process_gene_gene_from_web']

import pandas
import logging
from .processor import GnbrGeneGeneProcessor

base_url = 'https://zenodo.org/record/3459420/files'
logger = logging.getLogger(__name__)


def process_gene_gene(part1_path, part2_path):
    logger.info(f'Loading part 1 table from {part1_path}')
    df1 = pandas.read_csv(part1_path, sep='\t')
    logger.info(f'Loading part 2 table from {part2_path}')
    df2 = pandas.read_csv(part2_path, sep='\t', header=None)
    gp = GnbrGeneGeneProcessor(df1, df2)
    gp.extract_activations()
    return gp


def process_gene_gene_from_web():
    fname1 = f'{base_url}/part-i-gene-gene-path-theme-distributions.txt.gz'
    fname2 = f'{base_url}/part-ii-dependency-paths-gene-gene-sorted-with-themes.txt.gz'
    return process_gene_gene(fname1, fname2)