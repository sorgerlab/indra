__all__ = ['process_gene_gene', 'process_gene_gene_from_web',
           'process_gene_disease', 'process_gene_disease_from_web',
           'process_chemical_disease', 'process_chemical_disease_from_web',
           'process_chemical_gene', 'process_chemical_gene_from_web']

import pandas as pd
import logging
from .processor import GnbrProcessor

base_url = 'https://zenodo.org/record/3459420/files'
logger = logging.getLogger(__name__)


def process_gene_gene(part1_path: str, part2_path: str) -> GnbrProcessor:
    """Process gene–gene interactions.

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
        A GnbrProcessor object which contains a list of extracted INDRA
        Statements in its statements attribute.
    """
    return _process_from_files(part1_path, part2_path, 'gene', 'gene')


def process_chemical_gene(part1_path: str, part2_path: str) -> GnbrProcessor:
    """Process chemical–gene interactions.

    Parameters
    ----------
    part1_path :
        Path to the first dataset of dependency paths and themes.
    part2_path :
        Path to the second dataset of dependency paths and entity pairs.

    Returns
    -------
    gp :
        A GnbrProcessor object which contains a list of extracted
        INDRA Statements in its statements attribute.
    """
    return _process_from_files(part1_path, part2_path, 'chemical', 'gene')


def process_gene_disease(part1_path: str, part2_path: str) -> GnbrProcessor:
    """Process gene–disease interactions.

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
        A GnbrProcessor object which contains a list of extracted INDRA
        Statements in its statements attribute.
    """
    return _process_from_files(part1_path, part2_path, 'gene', 'disease')


def process_chemical_disease(part1_path: str, part2_path: str) \
        -> GnbrProcessor:
    """Process chemical–disease interactions.

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
        A GnbrProcessor object which contains a list of extracted INDRA
        Statements in its statements attribute.
    """
    return _process_from_files(part1_path, part2_path, 'chemical', 'disease')


def _process_from_files(part1_path, part2_path, first_type, second_type):
    logger.info(f'Loading part 1 table from {part1_path}')
    df1: pd.DataFrame = pd.read_csv(part1_path, sep='\t')
    logger.info(f'Loading part 2 table from {part2_path}')
    df2: pd.DataFrame = pd.read_csv(part2_path, sep='\t', header=None)
    gp: GnbrProcessor = GnbrProcessor(df1, df2, first_type, second_type)
    gp.extract_stmts()
    return gp


def process_gene_gene_from_web() -> GnbrProcessor:
    """Call process_gene_gene function on the GNBR datasets."""
    return _process_from_web('gene', 'gene')


def process_chemical_gene_from_web() -> GnbrProcessor:
    """Call process_chemical_gene function on the GNBR datasets."""
    return _process_from_web('chemical', 'gene')


def process_gene_disease_from_web() -> GnbrProcessor:
    """Call process_gene_disease function on the GNBR datasets."""
    return _process_from_web('gene', 'disease')


def process_chemical_disease_from_web() -> GnbrProcessor:
    """Call process_chemical_disease function on the GNBR datasets."""
    return _process_from_web('chemical', 'disease')


def _process_from_web(first_type, second_type):
    fname1 = (f'{base_url}/part-i-{first_type}-{second_type}-path-theme-'
              f'distributions.txt.gz')
    fname2 = (f'{base_url}/part-ii-dependency-paths-{first_type}-{second_type}'
              f'-sorted-with-themes.txt.gz')
    return _process_from_files(fname1, fname2, first_type, second_type)