import requests
import pandas as pd
import logging
# Python3
try:
    from io import StringIO
# Python2
except ImportError:
    from StringIO import StringIO

logger = logging.getLogger('cbio')
cbio_url = 'http://www.cbioportal.org/webservice.do'

ccle_study = 'cellline_ccle_broad'


def send_request(data, skiprows=0):
    '''
    Sends a web service requrest to the cBio portal with arguments given in
    the dictionary data and returns a Pandas data frame on success.

    Parameters
    ----------
    data : dict
        A dict of parameters for the cBioPortal query
    skiprows : int
        Number of rows to skip when reading dataframe. This is useful to align
        headers

    Returns
    -------
    df : Pandas DataFrame
        return the response from cBioPortal as a Pandas DataFrame
    '''
    res = requests.get(cbio_url, params=data)
    status = res.status_code
    if status == 200:
        csv_StringIO = StringIO(res.text)
        df = pd.read_csv(csv_StringIO, sep='\t', skiprows=skiprows)
        return df
    else:
        logger.error('Request returned with code %d' % status)


def get_mutations_ccle_lines_genes(lines, gene_list):
    '''
    Given a list of cell lines and genes, return the mutations in
    those lines and genes, if any

    Parameters
    ----------
    lines : list of str
        the names of the CCLE cell line(s)
    gene_list : list of str
        HGNC names of gene(s)

    Returns
    -------
    mutations : dict
        return the response from cBioPortal as a dict in the format
        {cell_line : {gene : [mutation1, mutation2, ...] }}

        for example -
        {'LOXIMVI': {'BRAF': ['V600E', 'I208V']},
         'SMEL30': {'BRAF': ['V600E', 'I208V']}}
    '''
    gene_str = ''
    for x in gene_list:
        gene_str += (x + ', ')
    gene_str[:len(gene_str) - 2]
    data = {'cmd': 'getMutationData',
            'case_set_id': ccle_study,
            'genetic_profile_id': ccle_study + '_mutations',
            'gene_list': gene_str}
    df = send_request(data, skiprows=1)
    df_lines_dict = {x.split('_')[0]: x
                     for x in df['case_id'].unique().tolist()}
    filter_lines = [df_lines_dict.get(x, None) for x in lines]
    filter_lines = [x for x in filter_lines if x is not None]
    df = df[df['case_id'].isin(filter_lines)]
    mutations = {}
    for c in lines:
        line_mutations = {}
        for g in gene_list:
            df_g = df[df['gene_symbol'] == g]
            amino_acid_change = df_g['amino_acid_change'].tolist()
            line_mutations[g] = amino_acid_change
        mutations[c] = line_mutations
    return mutations


def check_ccle_lines_for_mutation(gene, amino_acid_change):
    '''
    Check which cell lines in CCLE have a particular mutation

    Parameters
    ----------
    gene: str
        as HGNC ID
    amino_acid_change: str
        example - V600E

    Returns
    -------
    cell_lines : list
        return the response from cBioPortal as a list of cell lines
    '''
    data = {'cmd': 'getMutationData',
            'case_set_id': ccle_study,
            'genetic_profile_id': ccle_study + '_mutations',
            'gene_list': gene}
    df = send_request(data, skiprows=1)
    df = df[df['amino_acid_change'] == amino_acid_change]
    cell_lines = df['case_id'].unique().tolist()
    cell_lines = [x.split('_')[0] for x in cell_lines]
    return cell_lines
