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


def process_from_web(subset, url=None):
    if subset not in urls:
        raise ValueError('%s is not a valid CTD subset.' % subset)
    url = url if url else urls[subset]
    return _process_url_or_file(url, subset)


def process_tsv(fname, subset):
    return _process_url_or_file(fname, subset)


def _process_url_or_file(path, subset):
    df = pandas.read_csv(path, sep='\t', comment='#',
                         header=None, dtype=str, keep_default_na=False)
    return process_dataframe(df, subset)


def process_dataframe(df, subset):
    if subset not in processors:
        raise ValueError('%s is not a valid CTD subset.' % subset)
    cp = processors[subset](df)
    cp.extract_statements()
    return cp
