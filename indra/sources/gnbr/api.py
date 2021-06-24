__all__ = ['process_gene_gene', 'process_gene_gene_from_web']

import pandas as pd
import logging
from .processor import GnbrGeneGeneProcessor

base_url = 'https://zenodo.org/record/3459420/files'
logger = logging.getLogger(__name__)


def process_gene_gene(part1_path: str, part2_path: str) \
        -> GnbrGeneGeneProcessor:
    """Process gene-gene interactions.

    Parameters
    ----------
    part1_path :
        Path to the first dataset which contains dependency paths and themes.
    part2_path :
        Path to the second dataset which contains dependency paths and entity
        pairs.

    Returns
    -------
    gp :
        A GnbrGeneGeneProcessor object which contains a list of extracted INDRA
        Statements in its statements attribute.
    """
    logger.info(f'Loading part 1 table from {part1_path}')
    df1: pd.DataFrame = pd.read_csv(part1_path, sep='\t')
    logger.info(f'Loading part 2 table from {part2_path}')
    df2: pd.DataFrame = pd.read_csv(part2_path, sep='\t', header=None)
    gp: GnbrGeneGeneProcessor = GnbrGeneGeneProcessor(df1, df2)
    gp.extract_activations()
    return gp


def process_gene_gene_from_web() -> GnbrGeneGeneProcessor:
    """
    Call the process_gene_gene function with the paths to the GNBR datasets.
    """
    fname1 = f'{base_url}/part-i-gene-gene-path-theme-distributions.txt.gz'
    fname2 = (f'{base_url}/part-ii-dependency-paths-gene-gene-sorted-with-'
              'themes.txt.gz')
    return process_gene_gene(fname1, fname2)
