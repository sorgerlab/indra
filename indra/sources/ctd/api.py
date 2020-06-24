import pandas
from .processor import CTDProcessor


def process_from_web():
    url = 'http://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz'
    df = pandas.read_csv(url, sep='\t', skiprows=29)
    return process_dataframe(df)


def process_tsv(fname):
    df = pandas.read_csv(fname, sep='\t', skiprows=29)
    return process_dataframe(df)


def process_dataframe(df):
    cp = CTDProcessor(df)
    cp.extract_chemical_gene()
    return cp
