__all__ = ['process_from_web', 'process_files', 'process_df']


import pandas
from .processor import AcsnProcessor

ACSN_URL = 'https://acsn.curie.fr/ACSN2/downloads/'
ACSN_RELATIONS_URL = ACSN_URL + \
    'ACSN2_binary_relations_between_proteins_with_PMID.txt'
ACSN_CORRESPONDENCE_URL = ACSN_URL + 'ACSN2_HUGO_Correspondence.gmt'


def process_from_web():
    relations_df = pandas.read_csv(ACSN_RELATIONS_URL)
    correspondence_df = pandas.read_csv(ACSN_CORRESPONDENCE_URL)
    return process_df(relations_df, correspondence_df)


def process_files(relations_path: str, correspondence_path: str):
    relations_df = pandas.read_csv(relations_path)
    correspondence_df = pandas.read_csv(correspondence_path)
    return process_df(relations_df, correspondence_df)


def process_df(relations_df, correspondence_df):
    ap = AcsnProcessor(relations_df, correspondence_df)
    ap.extract_statements()
    return ap
