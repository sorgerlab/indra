__all__ = ['process_from_web', 'process_files', 'process_df']

from typing import Mapping

import pandas
import requests
from .processor import AcsnProcessor

ACSN_URL = 'https://acsn.curie.fr/ACSN2/downloads/'
ACSN_RELATIONS_URL = ACSN_URL + \
                     'ACSN2_binary_relations_between_proteins_with_PMID.txt'
ACSN_CORRESPONDENCE_URL = ACSN_URL + 'ACSN2_HUGO_Correspondence.gmt'


def process_from_web() -> AcsnProcessor:
    """Process ACSN data directly from the web.

    Returns
    -------
    :
        A processor with a list of INDRA statements that were extracted
        in its statements attribute.
    """
    relations_df = pandas.read_csv(ACSN_RELATIONS_URL, sep='\t')
    correspondence_dict = _transform_gmt(
        requests.get(ACSN_CORRESPONDENCE_URL).text.split('\n'))
    return process_df(relations_df, correspondence_dict)


def process_files(relations_path: str, correspondence_path: str) -> \
        AcsnProcessor:
    """Process ACSN data from input files.

    Parameters
    ----------
    relations_path :
        Path to the ACSN binary relations file.
    correspondence_path :
        Path to the ACSN correspondence GMT file.

    Returns
    -------
    :
        A processor with a list of INDRA statements that were extracted
        in its statements attribute.
    """
    relations_df = pandas.read_csv(relations_path)
    with open(correspondence_path, 'r') as fh:
        correspondence_dict = _transform_gmt(fh)
    return process_df(relations_df, correspondence_dict)


def process_df(relations_df: pandas.DataFrame, correspondence_dict: Mapping) \
        -> AcsnProcessor:
    """Process ACSN data from input data structures.

    Parameters
    ----------
    relations_df :
        An ACSN tab-separated data frame which consists of binary relationships
        between proteins with PMIDs.
    correspondence_dict :
        A dictionary with correspondences between ACSN entities and their HGNC
        symbols.

    Returns
    -------
    :
        A processor with a list of INDRA statements that were extracted
        in its statements attribute.
    """
    ap = AcsnProcessor(relations_df, correspondence_dict)
    ap.extract_statements()
    return ap


def _transform_gmt(gmt):
    """Convert ACSN correspondence GMT file into a dictionary."""
    acsn_hgnc_dict = {}
    for line in gmt:
        parts = line.strip().split('\t')
        acsn_hgnc_dict[parts[0]] = parts[2:]
    return acsn_hgnc_dict
