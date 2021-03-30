# -*- coding: utf-8 -*-

"""Code for downloading `Drug Gene Interaction Database (DGI-DB) <http://www.dgidb.org>`_."""

import logging
from typing import Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

USECOLS = [
    'gene_name', 'entrez_id', 'interaction_claim_source',
    'interaction_types', 'drug_name', 'drug_concept_id', 'PMIDs',
]


def get_version_df(version: Optional[str] = None) -> Tuple[str, pd.DataFrame]:
    """Get the latest version of the DGI interaction dataframe."""
    if version is None:
        try:
            import bioversions
            version = bioversions.get_version('Drug Gene Interaction Database')
        except ImportError:
            version = None
    if version is None:
        logger.warning('could not find version with bioregistry')
        version = '2021-Jan'
    url = f'https://www.dgidb.org/data/monthly_tsvs/{version}/interactions.tsv'
    df = pd.read_csv(url, usecols=USECOLS, sep='\t', dtype=str)
    return version, df
