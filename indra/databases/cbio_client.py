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
    """Return a data frame from a web service request to cBio portal.

    Sends a web service requrest to the cBio portal with arguments given in
    the dictionary data and returns a Pandas data frame on success.

    More information about the service here:
    http://www.cbioportal.org/web_api.jsp

    Parameters
    ----------
    data : dict
        A dict of parameters for the query.
    skiprows : int
        Number of rows to skip when reading dataframe. This is useful to align
        headers.

    Returns
    -------
    df : Pandas DataFrame
        return the response from cBioPortal as a Pandas DataFrame
    """
    res = requests.get(cbio_url, params=data)
    if res.status_code == 200:
        csv_StringIO = StringIO(res.text)
        df = pd.read_csv(csv_StringIO, sep='\t', skiprows=skiprows)
        return df
    else:
        logger.error('Request returned with code %d' % status)


def get_mutations_ccle_lines_genes(lines, gene_list):
    """Return a dict of mutations in given genes and cell lines.

    Parameters
    ----------
    lines : list[str]
        A list of CCLE cell line names to get mutations for
    gene_list : list[str]
        A list of HGNC gene symbols to get mutations in

    Returns
    -------
    mutations : dict
        The result from cBioPortal as a dict in the format
        {cell_line : {gene : [mutation1, mutation2, ...] }}

        Example:
        {'LOXIMVI': {'BRAF': ['V600E', 'I208V']},
         'SMEL30': {'BRAF': ['V600E', 'I208V']}}
    """
    gene_str = ','.join(gene_list)
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
        df_c = df[df['case_id'] == df_lines_dict[c]]
        line_mutations = {}
        for g in gene_list:
            df_g = df_c[df_c['gene_symbol'] == g]
            amino_acid_change = df_g['amino_acid_change'].tolist()
            line_mutations[g] = amino_acid_change
        mutations[c] = line_mutations
    return mutations


def check_ccle_lines_for_mutation(gene, amino_acid_change):
    """Return cell lines with a given mutation.

    Check which cell lines in CCLE have a particular mutation and
    return their names in a list.

    Parameters
    ----------
    gene : str
        as HGNC ID
    amino_acid_change : str
        example - V600E

    Returns
    -------
    cell_lines : list
        return the response from cBioPortal as a list of cell lines
    """
    data = {'cmd': 'getMutationData',
            'case_set_id': ccle_study,
            'genetic_profile_id': ccle_study + '_mutations',
            'gene_list': gene}
    df = send_request(data, skiprows=1)
    df = df[df['amino_acid_change'] == amino_acid_change]
    cell_lines = df['case_id'].unique().tolist()
    cell_lines = [x.split('_')[0] for x in cell_lines]
    return cell_lines
