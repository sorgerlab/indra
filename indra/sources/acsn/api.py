__all__ = ['process_from_web', 'process_df']

import pandas
import requests
from .processor import AcsnProcessor

ACSN_URL = 'https://acsn.curie.fr/ACSN2/downloads/'
ACSN_RELATIONS_URL = ACSN_URL + \
                     'ACSN2_binary_relations_between_proteins_with_PMID.txt'
ACSN_CORRESPONDENCE_URL = ACSN_URL + 'ACSN2_HUGO_Correspondence.gmt'


def process_from_web():
    """Process ACSN relations and correspondence data files
    from web into INDRA statements.

    Returns
    -------
    AcsnProcessor
        A ACSNProcessor which contains INDRA statements extracted from
        given ACSN relations and correspondence gmt
    """
    relations_df = pandas.read_csv(ACSN_RELATIONS_URL, sep='\t')
    correspondence_dict = _transform_gmt(
        requests.get(ACSN_CORRESPONDENCE_URL).text.split('\n'))
    return process_df(relations_df, correspondence_dict)


def process_files(relations_path: str, correspondence_path: str):
    relations_df = pandas.read_csv(relations_path)
    with open(correspondence_path, 'r') as fh:
        correspondence_dict = _transform_gmt(fh)
    return process_df(relations_df, correspondence_dict)


def process_df(relations_df, correspondence_dict):
    """Get ACSNProcessor which extracted INDRA statements from given
    ACSN relations dataframe

    Parameters
    ----------
    relations_df : pandas.DataFrame
        A tab-separated data frame which consists of binary relationship between
        proteins with PMIDs
    correspondence_dict : dict
        A dictionary with Correspondence between ACSN entities and its HUGO names

    Returns
    -------
    ap : AcsnProcessor
        A processor with a list of INDRA statements that were extracted
    """
    ap = AcsnProcessor(relations_df, correspondence_dict)
    ap.extract_statements()
    return ap


def _transform_gmt(gmt):
    """Convert ACSN correspondence GMT file into a dictionary"""
    acsn_hgnc_dict = {}
    for line in gmt:
        parts = line.strip().split('\t')
        acsn_hgnc_dict[parts[0]] = parts[2:]
    return acsn_hgnc_dict
