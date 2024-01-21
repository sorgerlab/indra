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

import logging
import requests
from functools import lru_cache
from indra.databases import hgnc_client


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
    method : str
        The HTTP method to use for the request.
        Example: 'get' or 'post'
    endpoint : str
        The endpoint to use for the request.
        Example: 'studies'
    json_data : Optional[Dict]
        The dict-like JSON data structure to send with the request.

    Returns
    -------
    JSON
        The JSON object returned by the web service call.
    """
    if endpoint.startswith('/'):
        endpoint = endpoint[1:]
    request_fun = getattr(requests, method)
    full_url = cbio_url + '/' + endpoint
    print('URL: %s' % full_url)
    print('JSON: %s' % json_data)
    res = request_fun(full_url, json=json_data or {})
    if res.status_code != 200:
        logger.error(f'Request returned with code {res.status_code}: '
                     f'{res.text}')
        return
    return res.json()


def get_mutations(study_id, gene_list=None, mutation_type=None,
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
    mutations : dict
        A dict with entries for each gene symbol and another list
        with entries for each corresponding amino acid change.
    """
    genetic_profile = get_genetic_profiles(study_id, 'mutation')[0]

    entrez_to_gene_symbol = get_entrez_mappings(gene_list)
    entrez_ids = list(entrez_to_gene_symbol)

    # Does this need to be parameterized?
    case_set_id = study_id + '_all'

    mutations = send_request('post',
                             f'molecular-profiles/{genetic_profile}/'
                             f'mutations/fetch',
                             {'sampleListId': case_set_id,
                              'entrezGeneIds': entrez_ids})

    if case_id:
        mutations = [m for m in mutations if m['sampleId'] == case_id]

    if mutation_type:
        mutations = [m for m in mutations if (mutation_type.casefold()
                                              in m['mutationType'].casefold())]

    mutations_dict = {
        'gene_symbol': [entrez_to_gene_symbol[str(m['entrezGeneId'])]
                        for m in mutations],
        'amino_acid_change': [m['proteinChange'] for m in mutations],
        'sample_id': [m['sampleId'] for m in mutations],
    }
    return mutations_dict


def get_entrez_mappings(gene_list):
    if gene_list:
        # First we need to get HGNC IDs from HGNC symbols
        hgnc_mappings = {g: hgnc_client.get_hgnc_id(g) for g in gene_list}
        # Next, we map from HGNC symbols to Entrez IDs via the hgnc_mappings
        entrez_mappings = {g: hgnc_client.get_entrez_id(hgnc_mappings[g])
                           for g in gene_list if hgnc_mappings[g] is not None}
        # Finally, we reverse the mapping, this will ensure that
        # we can get the gene symbols back when generating results
        entrez_to_gene_symbol = {v: k for k, v in entrez_mappings.items()
                                 if v is not None and k is not None}
    else:
        entrez_to_gene_symbol = {}
    return entrez_to_gene_symbol


def get_case_lists(study_id):
    """Return a list of the case set ids for a particular study.

    In v2 of the API these are called sample lists.

    Old comment:
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
    case_set_ids : list[str]
        A list of case set IDs, e.g., ['cellline_ccle_broad_all',
        'cellline_ccle_broad_cna', ...]
    """
    res = send_request('get', f'studies/{study_id}/sample-lists')
    return [sl['sampleListId'] for sl in res]


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
        A string that specifies which case_set_id to use, based on a complete
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
    case_set_ids = get_case_lists(study_id)
    if case_set_filter:
        case_set_id = [x for x in case_set_ids if case_set_filter in x][0]
    else:
        # based on looking at the cBioPortal, this is a common case_set_id
        case_set_id = study_id + '_all'
    entrez_to_gene_symbol = get_entrez_mappings(gene_list)
    entrez_ids = list(entrez_to_gene_symbol)
    res = send_request('post', f'molecular-profiles/{genetic_profile}/'
                               f'molecular-data/fetch',
                       {'sampleListId': case_set_id,
                        'entrezGeneIds': entrez_ids})

    profile_data = {}
    # Each entry in the results contains something like
    # {'entrezGeneId': 673, 'molecularProfileId': 'cellline_ccle_broad_cna',
    #  'sampleId': '1321N1_CENTRAL_NERVOUS_SYSTEM',
    #  'studyId': 'cellline_ccle_broad', 'value': 1, ...}
    for sample in res:
        sample_id = sample['sampleId']
        if sample_id not in profile_data:
            profile_data[sample_id] = {}
        gene_symbol = entrez_to_gene_symbol[str(sample['entrezGeneId'])]
        profile_data[sample_id][gene_symbol] = sample['value']
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
    # First we get all the case lists for the study
    case_lists = get_case_lists(study_id)
    # Then we find ones that have 'sequenced' in the name
    sequencing_case_list = [cl for cl in case_lists if 'sequenced' in cl]
    # Then we look at the sample IDs and count them
    cases = set()
    for cl in sequencing_case_list:
        res = send_request('get', f'/sample-lists/{cl}/sample-ids')
        cases |= set(res)
    num_case = len(cases)
    return num_case


def get_genetic_profiles(study_id, profile_filter=None):
    """Return all the genetic profiles (data sets) for a given study.

    Genetic profiles are different types of data for a given study. For
    instance the study 'cellline_ccle_broad' has profiles such as
    'cellline_ccle_broad_mutations' for mutations, 'cellline_ccle_broad_CNA'
    for copy number alterations, etc.

    NOTE: In the v2 API, the genetic profiles are called molecular profiles.

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
    cancer_types = send_request('get', 'cancer-types')
    if cancer_filter:
        cancer_types = [c for c in cancer_types
                        if cancer_filter.casefold() in c['name'].casefold()]
    type_ids = [c['cancerTypeId'] for c in cancer_types]
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
    mutations = get_mutations(ccle_study, [gene], 'missense')
    cell_lines = {cl for aac, cl
                  in zip(mutations['amino_acid_change'], mutations['sample_id'])
                  if aac == amino_acid_change}
    return sorted(cell_lines)


def get_ccle_cna(gene_list, cell_lines=None):
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
    cell_lines : Optional[list[str]]
        A list of CCLE cell line names to get mutations for.

    Returns
    -------
    profile_data : dict[dict[int]]
        A dict keyed to cases containing a dict keyed to genes
        containing int
    """
    profile_data = get_profile_data(ccle_study, gene_list,
                                    'COPY_NUMBER_ALTERATION', 'all')
    return {cell_line: value for cell_line, value in profile_data.items()
            if cell_lines is None or cell_line in cell_lines}


def get_ccle_mrna(gene_list, cell_lines=None):
    """Return a dict of mRNA amounts in given genes and cell lines from CCLE.

    Parameters
    ----------
    gene_list : list[str]
        A list of HGNC gene symbols to get mRNA amounts for.
    cell_lines : Optional[list[str]]
        A list of CCLE cell line names to get mRNA amounts for.

    Returns
    -------
    mrna_amounts : dict[dict[float]]
        A dict keyed to cell lines containing a dict keyed to genes
        containing float
    """
    profile_data = get_profile_data(ccle_study, gene_list,
                                    'MRNA_EXPRESSION', 'all')
    # FIXME: we need a data structure like this
    #         assert mrna['A375_SKIN'] is not None
    #         assert mrna['A375_SKIN']['MAP2K1'] > 10
    # >       assert mrna['A375_SKIN']['XYZ'] is None
    # E       KeyError: 'XYZ'
    mrna_amounts = {cell_line: value
                    for cell_line, value in profile_data.items()
                    if cell_lines is None or cell_line in cell_lines}
    # This is to make sure that if cell_lines were specified then
    # we return None if there is no data for a given cell line
    # This matches the old behavior of the function
    if cell_lines:
        for cell_line in cell_lines:
            if cell_line not in mrna_amounts:
                mrna_amounts[cell_line] = None
    return mrna_amounts
