# -*- coding: utf-8 -*-

"""API for `Drug Gene Interaction DB <http://www.dgidb.org>`_."""

import logging
from typing import Optional, Tuple

import pandas as pd

from .processor import DGIProcessor

logger = logging.getLogger(__name__)

USECOLS = [
    "gene_name",
    "entrez_id",
    "interaction_claim_source",
    "interaction_types",
    "drug_name",
    "drug_concept_id",
    "PMIDs",
]


def process_version(version: Optional[str] = None) -> DGIProcessor:
    """Get a processor that extracted INDRA Statements from DGI content.

    Parameters
    ----------
    version :
        The optional version of DGI to use. If no ``df`` is given, this is
        also automatically looked up.

    Returns
    -------
    dp :
        A DGI processor with pre-extracted INDRA statements
    """
    version, df = get_version_df(version)
    return process_df(df=df, version=version)


def process_df(
    df: pd.DataFrame,
    version: Optional[str] = None,
) -> DGIProcessor:
    """Get a processor that extracted INDRA Statements from DGI content based
    on the given dataframe.

    Parameters
    ----------
    df :
        A pandas DataFrame for the DGI interactions file.
    version :
        The optional version of DGI to use. If not given, statements will
        not be annotated with a version number.

    Returns
    -------
    dp :
        A DGI processor with pre-extracted INDRA statements
    """
    dp = DGIProcessor(df=df, version=version)
    dp.extract_statements()
    return dp


def get_version_df(version: Optional[str] = None) -> Tuple[str, pd.DataFrame]:
    """Get the latest version of the DGI interaction dataframe."""
    if version is None:
        try:
            import bioversions
        except ImportError:
            version = None
        else:
            version = bioversions.get_version("Drug Gene Interaction Database")
    if version is None:
        logger.warning("could not find version with bioregistry")
        version = "2021-Jan"
    url = f"https://www.dgidb.org/data/monthly_tsvs/{version}/interactions.tsv"
    df = pd.read_csv(url, usecols=USECOLS, sep="\t", dtype=str)
    return version, df
