# -*- coding: utf-8 -*-

"""Code for downloading `Drug Gene Interaction Database (DGI-DB) <http://www.dgidb.org>`_."""

from typing import Optional, Tuple

import pandas as pd

USECOLS = [
    'gene_name', 'entrez_id', 'interaction_claim_source',
    'interaction_types', 'drug_name', 'drug_concept_id', 'PMIDs',
]


def get_version_df(version: Optional[str] = None) -> Tuple[str, pd.DataFrame]:
    """Get the latest version of the DGI interaction dataframe."""
    if version is None:
        version = '2021-Jan'  # TODO use bioversions to look up with following two lines
        # import bioregistry
        # version = bioregistry.get_version('dgi')
    url = f'https://www.dgidb.org/data/monthly_tsvs/{version}/interactions.tsv'
    df = pd.read_csv(url, usecols=USECOLS, dtype=str)
    return version, df
