"""This is a client for the cBioPortal web service, with
documentation at https://docs.cbioportal.org/web-api-and-clients/
and Swagger definition at https://www.cbioportal.org/api/v2/api-docs.
Note that the client implements direct requests to the API instead of
adding an additional dependency to do so.
"""
__all__ = ["get_mutations", "get_case_lists", "get_profile_data",
           "get_num_sequenced", "get_genetic_profiles",
           "get_cancer_studies", "get_cancer_types", "get_ccle_mutations",
           "get_ccle_lines_for_mutation", "get_ccle_cna",
           "get_ccle_mrna"]

import pandas
import logging
import requests


logger = logging.getLogger(__name__)

cbio_url = 'https://www.cbioportal.org/api'
ccle_study = 'cellline_ccle_broad'

# TODO: implement caching with json_data made immutable
def send_request(method, endpoint, json_data=None):
    """Return the results of a web service request to cBio portal.

    Sends a web service request to the cBio portal with a specific endpoint,
    method, and JSON data structure, and returns the resulting JSON
    data structure on success.

    More information about the service is available here:
    https://www.cbioportal.org/api/v2/api-docs

    Parameters
    ----------
    TODO

    Returns
    -------
    JSON
        The JSON object returned by the web service call.
    """
    if endpoint.startswith('/'):
        endpoint = endpoint[1:]
    request_fun = getattr(requests, method)
    res = request_fun(cbio_url + '/' + endpoint, json=json_data)
    if res.status_code != 200:
        logger.error('Request returned with code %d' % res.status_code)
        return
    return res.json()


def get_mutations(study_id, gene_list, mutation_type=None,
                  case_id=None):
    """Return mutations as a list of genes and list of amino acid changes.

    Parameters
    ----------
    study_id : str
        The ID of the cBio study.
        Example: 'cellline_ccle_broad' or 'paad_icgc'
    gene_list : list[str]
        A list of genes with their HGNC symbols.
        Example: ['BRAF', 'KRAS']
    mutation_type : Optional[str]
        The type of mutation to filter to.
        mutation_type can be one of: missense, nonsense, frame_shift_ins,
        frame_shift_del, splice_site
    case_id : Optional[str]
        The case ID within the study to filter to.

    Returns
    -------
    mutations : tuple[list]
        A tuple of two lists, the first one containing a list of genes, and
        the second one a list of amino acid changes in those genes.
    """
    genetic_profile = get_genetic_profiles(study_id, 'mutation')[0]
    gene_list_str = ','.join(gene_list)

    data = {'cmd': 'getMutationData',
            'case_set_id': study_id,
            'genetic_profile_id': genetic_profile,
            'gene_list': gene_list_str,
            'skiprows': -1}
    df = send_request(**data)
    if case_id:
        df = df[df['case_id'] == case_id]
    res = _filter_data_frame(df, ['gene_symbol', 'amino_acid_change'],
                             'mutation_type', mutation_type)
    mutations = {'gene_symbol': list(res['gene_symbol'].values()),
                 'amino_acid_change': list(res['amino_acid_change'].values())}
    return mutations


def get_case_lists(study_id):
    """Return a list of the case set ids for a particular study.

    TAKE NOTE the "case_list_id" are the same thing as "case_set_id"
    Within the data, this string is referred to as a "case_list_id".
    Within API calls it is referred to as a 'case_set_id'.
    The documentation does not make this explicitly clear.

    Parameters
    ----------
    study_id : str
        The ID of the cBio study.
        Example: 'cellline_ccle_broad' or 'paad_icgc'

    Returns
    -------
    case_set_ids : dict[dict[int]]
        A dict keyed to cases containing a dict keyed to genes
        containing int
    """
    res = send_request('post', 'sample-lists/fetch',
                       {'studyIds': [study_id]})
    data = {'cmd': 'getCaseLists',
            'cancer_study_id': study_id}
    df = send_request(**data)
    case_set_ids = df['case_list_id'].tolist()
    return case_set_ids


def get_profile_data(study_id, gene_list,
                     profile_filter, case_set_filter=None):
    """Return dict of cases and genes and their respective values.

    Parameters
    ----------
    study_id : str
        The ID of the cBio study.
        Example: 'cellline_ccle_broad' or 'paad_icgc'
    gene_list : list[str]
        A list of genes with their HGNC symbols.
        Example: ['BRAF', 'KRAS']
    profile_filter : str
        A string used to filter the profiles to return. Will be one of:
        - MUTATION
        - MUTATION_EXTENDED
        - COPY_NUMBER_ALTERATION
        - MRNA_EXPRESSION
        - METHYLATION
    case_set_filter : Optional[str]
        A string that specifices which case_set_id to use, based on a complete
        or partial match. If not provided, will look for study_id + '_all'

    Returns
    -------
    profile_data : dict[dict[int]]
        A dict keyed to cases containing a dict keyed to genes
        containing int
    """
    genetic_profiles = get_genetic_profiles(study_id, profile_filter)
    if genetic_profiles:
        genetic_profile = genetic_profiles[0]
    else:
        return {}
    gene_list_str = ','.join(gene_list)
    case_set_ids = get_case_lists(study_id)
    if case_set_filter:
        case_set_id = [x for x in case_set_ids if case_set_filter in x][0]
    else:
        case_set_id = study_id + '_all'
        # based on looking at the cBioPortal, this is a common case_set_id
    data = {'cmd': 'getProfileData',
            'case_set_id': case_set_id,
            'genetic_profile_id': genetic_profile,
            'gene_list': gene_list_str,
            'skiprows': -1}
    df = send_request(**data)
    case_list_df = [x for x in df.columns.tolist()
                    if x not in ['GENE_ID', 'COMMON']]
    profile_data = {case: {g: None for g in gene_list}
                    for case in case_list_df}
    for case in case_list_df:
        profile_values = df[case].tolist()
        df_gene_list = df['COMMON'].tolist()
        for g, cv in zip(df_gene_list, profile_values):
            if not pandas.isnull(cv):
                profile_data[case][g] = cv
    return profile_data


def get_num_sequenced(study_id):
    """Return number of sequenced tumors for given study.

    This is useful for calculating mutation statistics in terms of the
    prevalence of certain mutations within a type of cancer.

    Parameters
    ----------
    study_id : str
        The ID of the cBio study.
        Example: 'paad_icgc'

    Returns
    -------
    num_case : int
        The number of sequenced tumors in the given study
    """
    data = {'cmd': 'getCaseLists',
            'cancer_study_id': study_id}
    df = send_request(**data)
    if df.empty:
        return 0
    row_filter = df['case_list_id'].str.contains('sequenced', case=False)
    num_case = len(df[row_filter]['case_ids'].tolist()[0].split(' '))
    return num_case


def get_genetic_profiles(study_id, profile_filter=None):
    """Return all the genetic profiles (data sets) for a given study.

    Genetic profiles are different types of data for a given study. For
    instance the study 'cellline_ccle_broad' has profiles such as
    'cellline_ccle_broad_mutations' for mutations, 'cellline_ccle_broad_CNA'
    for copy number alterations, etc.

    Parameters
    ----------
    study_id : str
        The ID of the cBio study.
        Example: 'paad_icgc'
    profile_filter : Optional[str]
        A string used to filter the profiles to return.
        Will be one of:
        - MUTATION
        - MUTATION_EXTENDED
        - COPY_NUMBER_ALTERATION
        - MRNA_EXPRESSION
        - METHYLATION
        The genetic profiles can include "mutation", "CNA", "rppa",
        "methylation", etc. The filter is case insensitive.

    Returns
    -------
    genetic_profiles : list[str]
        A list of genetic profiles available  for the given study.
    """
    res = send_request('get', f'studies/{study_id}/molecular-profiles')
    if profile_filter:
        res = [prof for prof in res
               if (profile_filter.casefold()
                   in prof['molecularAlterationType'].casefold())]
    profile_ids = [prof['molecularProfileId'] for prof in res]
    return profile_ids


def get_cancer_studies(study_filter=None):
    """Return a list of cancer study identifiers, optionally filtered.

    There are typically multiple studies for a given type of cancer and
    a filter can be used to constrain the returned list.

    Parameters
    ----------
    study_filter : Optional[str]
        A string used to filter the study IDs to return. Example: "paad"

    Returns
    -------
    study_ids : list[str]
        A list of study IDs.
        For instance "paad" as a filter would result in a list
        of study IDs with paad in their name like "paad_icgc", "paad_tcga",
        etc.
    """
    studies = send_request('get', 'studies')
    if study_filter:
        studies = [s for s in studies
                   if study_filter.casefold() in s['studyId'].casefold()]
    study_ids = [s['studyId'] for s in studies]
    return study_ids


def get_cancer_types(cancer_filter=None):
    """Return a list of cancer types, optionally filtered.

    Parameters
    ----------
    cancer_filter : Optional[str]
        A string used to filter cancer types. Its value is the name or
        part of the name of a type of cancer. Example: "melanoma",
        "pancreatic", "non-small cell lung"

    Returns
    -------
    type_ids : list[str]
        A list of cancer types matching the filter.
        Example: for cancer_filter="pancreatic", the result includes
        "panet" (neuro-endocrine) and "paad" (adenocarcinoma)
    """
    data = {'cmd': 'getTypesOfCancer'}
    df = send_request(**data)
    res = _filter_data_frame(df, ['type_of_cancer_id'], 'name', cancer_filter)
    type_ids = list(res['type_of_cancer_id'].values())
    return type_ids


def get_ccle_mutations(gene_list, cell_lines, mutation_type=None):
    """Return a dict of mutations in given genes and cell lines from CCLE.

    This is a specialized call to get_mutations tailored to CCLE cell lines.

    Parameters
    ----------
    gene_list : list[str]
        A list of HGNC gene symbols to get mutations in
    cell_lines : list[str]
        A list of CCLE cell line names to get mutations for.
    mutation_type : Optional[str]
        The type of mutation to filter to.
        mutation_type can be one of: missense, nonsense, frame_shift_ins,
        frame_shift_del, splice_site

    Returns
    -------
    mutations : dict
        The result from cBioPortal as a dict in the format
        {cell_line : {gene : [mutation1, mutation2, ...] }}

        Example:
        {'LOXIMVI_SKIN': {'BRAF': ['V600E', 'I208V']},
        'SKMEL30_SKIN': {'BRAF': ['D287H', 'E275K']}}
    """
    mutations = {cl: {g: [] for g in gene_list} for cl in cell_lines}
    for cell_line in cell_lines:
        mutations_cl = get_mutations(ccle_study, gene_list,
                                     mutation_type=mutation_type,
                                     case_id=cell_line)
        for gene, aa_change in zip(mutations_cl['gene_symbol'],
                                   mutations_cl['amino_acid_change']):
            aa_change = str(aa_change)
            mutations[cell_line][gene].append(aa_change)
    return mutations


def get_ccle_lines_for_mutation(gene, amino_acid_change):
    """Return cell lines with a given point mutation in a given gene.

    Checks which cell lines in CCLE have a particular point mutation
    in a given gene and return their names in a list.

    Parameters
    ----------
    gene : str
        The HGNC symbol of the mutated gene in whose product the amino
        acid change occurs. Example: "BRAF"
    amino_acid_change : str
        The amino acid change of interest. Example: "V600E"

    Returns
    -------
    cell_lines : list
        A list of CCLE cell lines in which the given mutation occurs.
    """
    data = {'cmd': 'getMutationData',
            'case_set_id': ccle_study,
            'genetic_profile_id': ccle_study + '_mutations',
            'gene_list': gene,
            'skiprows': 1}
    df = send_request(**data)
    df = df[df['amino_acid_change'] == amino_acid_change]
    cell_lines = df['case_id'].unique().tolist()
    return cell_lines


def get_ccle_cna(gene_list, cell_lines):
    """Return a dict of CNAs in given genes and cell lines from CCLE.

    CNA values correspond to the following alterations

    -2 = homozygous deletion

    -1 = hemizygous deletion

    0 = neutral / no change

    1 = gain

    2 = high level amplification

    Parameters
    ----------
    gene_list : list[str]
        A list of HGNC gene symbols to get mutations in
    cell_lines : list[str]
        A list of CCLE cell line names to get mutations for.

    Returns
    -------
    profile_data : dict[dict[int]]
        A dict keyed to cases containing a dict keyed to genes
        containing int
    """
    profile_data = get_profile_data(ccle_study, gene_list,
                                    'COPY_NUMBER_ALTERATION', 'all')
    profile_data = dict((key, value) for key, value in profile_data.items()
                        if key in cell_lines)
    return profile_data


def get_ccle_mrna(gene_list, cell_lines):
    """Return a dict of mRNA amounts in given genes and cell lines from CCLE.

    Parameters
    ----------
    gene_list : list[str]
        A list of HGNC gene symbols to get mRNA amounts for.
    cell_lines : list[str]
        A list of CCLE cell line names to get mRNA amounts for.

    Returns
    -------
    mrna_amounts : dict[dict[float]]
        A dict keyed to cell lines containing a dict keyed to genes
        containing float
    """
    gene_list_str = ','.join(gene_list)
    data = {'cmd': 'getProfileData',
            'case_set_id': ccle_study + '_mrna',
            'genetic_profile_id': ccle_study + '_mrna',
            'gene_list': gene_list_str,
            'skiprows': -1}
    df = send_request(**data)
    mrna_amounts = {cl: {g: [] for g in gene_list} for cl in cell_lines}
    for cell_line in cell_lines:
        if cell_line in df.columns:
            for gene in gene_list:
                value_cell = df[cell_line][df['COMMON'] == gene]
                if value_cell.empty:
                    mrna_amounts[cell_line][gene] = None
                elif pandas.isnull(value_cell.values[0]):
                    mrna_amounts[cell_line][gene] = None
                else:
                    value = value_cell.values[0]
                    mrna_amounts[cell_line][gene] = value
        else:
            mrna_amounts[cell_line] = None
    return mrna_amounts


def _filter_data_frame(df, data_col, filter_col, filter_str=None):
    """Return a filtered data frame as a dictionary."""
    if filter_str is not None:
        relevant_cols = data_col + [filter_col]
        df.dropna(inplace=True, subset=relevant_cols)
        row_filter = df[filter_col].str.contains(filter_str, case=False)
        data_list = df[row_filter][data_col].to_dict()
    else:
        data_list = df[data_col].to_dict()
    return data_list


# Deactivate this section for the time being, can be reinstated
# once these are fully integrated
'''

def _read_ccle_cna():
    fname = os.path.dirname(os.path.abspath(__file__)) + \
        '/../../data/ccle_CNA.txt'
    try:
        df = pandas.read_csv(fname, sep='\t')
    except Exception:
        df = None
    return df

ccle_cna_df = _read_ccle_cna()


def _read_ccle_mrna():
    fname = os.path.dirname(os.path.abspath(__file__)) + \
        '/../../data/ccle_expression_median.txt'
    try:
        df = pandas.read_csv(fname, sep='\t')
    except Exception:
        df = None
    return df

ccle_mrna_df = _read_ccle_mrna()


def _read_ccle_mutations():
    fname = os.path.dirname(os.path.abspath(__file__)) + \
        '/../../data/ccle_mutations_extended.txt'
    try:
        df = pandas.read_csv(fname, sep='\t', skiprows=2)
    except Exception:
        df = None
    return df

ccle_mutations_df = _read_ccle_mutations()
'''
