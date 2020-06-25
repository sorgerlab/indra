import pandas
from .processor import CTDProcessor, CTDChemicalDiseaseProcessor, \
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
    """Process a subset of CTD from the web into INDRA Statements.

    Parameters
    ----------
    subset : str
        A CTD subset, one of chemical_gene, chemical_disease,
        gene_disease.
    url : Optional[str]
        If not provided, the default CTD URL is used (beware, it usually
        gives permission denied). If provided, the given URL is used to
        access a tsv or tsv.gz file.

    Returns
    -------
    CTDProcessor
        A CTDProcessor which contains INDRA Statements extracted from the
        given CTD subset as its statements attribute.
    """
    if subset not in urls:
        raise ValueError('%s is not a valid CTD subset.' % subset)
    url = url if url else urls[subset]
    return _process_url_or_file(url, subset)


def process_tsv(fname, subset):
    """Process a subset of CTD from a tsv or tsv.gz file into INDRA Statements.

    Parameters
    ----------
    fname : str
        Path to a tsv or tsv.gz file of the given CTD subset.
    subset : str
        A CTD subset, one of chemical_gene, chemical_disease,
        gene_disease.

    Returns
    -------
    CTDProcessor
        A CTDProcessor which contains INDRA Statements extracted from the
        given CTD subset as its statements attribute.
    """
    return _process_url_or_file(fname, subset)


def _process_url_or_file(path, subset):
    df = pandas.read_csv(path, sep='\t', comment='#',
                         header=None, dtype=str, keep_default_na=False)
    return process_dataframe(df, subset)


def process_dataframe(df, subset):
    """Process a subset of CTD from a DataFrame into INDRA Statements.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame of the given CTD subset.
    subset : str
        A CTD subset, one of chemical_gene, chemical_disease,
        gene_disease.

    Returns
    -------
    CTDProcessor
        A CTDProcessor which contains INDRA Statements extracted from the
        given CTD subset as its statements attribute.
    """
    if subset not in processors:
        raise ValueError('%s is not a valid CTD subset.' % subset)
    cp = processors[subset](df)
    cp.extract_statements()
    return cp
