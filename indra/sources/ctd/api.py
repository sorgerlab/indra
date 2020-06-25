import pandas
from .processor import CTDChemicalDiseaseProcessor, \
    CTDGeneDiseaseProcessor, CTDChemicalGeneProcessor

base_url = 'http://ctdbase.org/reports/'

urls = {
    'chemical_gene': base_url + 'CTD_chem_gene_ixns.tsv.gz',
    'chemical_disease': base_url + 'CTD_chemicals_diseases.tsv.gz',
    'gene_disease': base_url + 'CTD_genes_diseases.tsv.gz',
}

processors = {
    'chemical_gene': CTDChemicalGeneProcessor,
    'chemical_disease': CTDChemicalDiseaseProcessor,
    'gene_disease': CTDGeneDiseaseProcessor,
}


def process_from_web(subset):
    if subset not in urls:
        raise ValueError('%s is not a valid CTD subset.')
    df = pandas.read_csv(urls[subset], sep='\t', comment='#',
                         header=None)
    return process_dataframe(df)


def process_tsv(fname, subset):
    df = pandas.read_csv(fname, sep='\t', comment='#', header=None)
    return process_dataframe(df, subset)


def process_dataframe(df, subset):
    if subset not in processors:
        raise ValueError('%s is not a valid CTD subset.')
    cp = processors[subset](df)
    cp.extract_statements()
    return cp
