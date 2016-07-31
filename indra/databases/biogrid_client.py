import os
import json
import logging
import requests

biogrid_url = 'http://webservice.thebiogrid.org/interactions/'

logger = logging.getLogger('biogrid')

# THIS FILE IS NOT UNDER VERSION CONTROL
# For more information see http://wiki.thebiogrid.org/doku.php/biogridrest
api_key_file = os.path.dirname(os.path.realpath(__file__)) + '/' + \
               'biogrid_api_key'
# Read the API key
try:
    with open(api_key_file, 'rt') as fh:
        api_key = fh.read().strip()
except IOError:
    logger.error('BioGRID API key could not be found.')
    logger.error(api_key_file)
    api_key = None

def get_publications(gene_names, save_json_name=None):
    """Return evidence publications for interaction between the given genes.

    Parameters
    ----------
    gene_names : list[str]
        A list of gene names (HGNC symbols) to query interactions between.
        Currently supports exactly two genes only.
    save_json_name : Optional[str]
        A file name to save the raw BioGRID web service output in. By default,
        the raw output is not saved.

    Return
    ------
    publications : list[Publication]
        A list of Publication objects that provide evidence for interactions
        between the given list of genes.
    """
    if len(gene_names) != 2:
        logger.warning('Other than 2 gene names given.')
        return []
    res_dict = _send_request(gene_names)
    if not res_dict:
        return []
    if save_json_name is not None:
        with open(save_json_name, 'wt') as fh:
            json.dump(res_dict, fh, indent=1)
    publications = _extract_publications(res_dict, gene_names)
    return publications

class Publication(object):
    def __init__(self, interaction, interaction_id):
        self.pmid = "PMID" + str(interaction['PUBMED_ID'])
        self.modification = interaction['MODIFICATION']
        self.experimental_system = interaction['EXPERIMENTAL_SYSTEM']
        self.experimental_system_type = interaction['EXPERIMENTAL_SYSTEM_TYPE']
        self.throughput = interaction['THROUGHPUT']
        self.interaction_id = interaction_id

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "Publication(%s)" % self.pmid

def _extract_publications(res_dict, gene_names):
    res_filtered = _filter_results(res_dict, gene_names)
    publications = []
    for interaction_id in res_filtered.keys():
        pub = Publication(res_filtered[interaction_id], interaction_id)
        publications.append(pub)
    return publications

def _filter_results(res_dict, gene_names):
    filtered_dict = {}
    for interaction_id in res_dict.keys():
        interactors = [res_dict[interaction_id]['OFFICIAL_SYMBOL_A'],
                       res_dict[interaction_id]['OFFICIAL_SYMBOL_B']]
        if set(interactors) == set(gene_names):
            filtered_dict[interaction_id] = res_dict[interaction_id]
    return filtered_dict

def _send_request(gene_names):
    if api_key is None:
        logger.error('BioGRID cannot be used without API key')
        return None
    params = {'searchNames': 'true',
              'geneList': '|'.join(gene_names),
              'taxId': '9606',
              'format': 'json',
              'includeInteractors': False,
              'accesskey': api_key}
    res = requests.get(biogrid_url, params)
    res.raise_for_status()
    res_dict = res.json()
    return res_dict
