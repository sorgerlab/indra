from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import os
import json
import logging
import requests
from collections import Counter
from indra.statements import Complex, Agent, Evidence

biogrid_url = 'http://webservice.thebiogrid.org/interactions/'

logger = logging.getLogger('biogrid')

# THIS FILE IS NOT UNDER VERSION CONTROL
# For more information see http://wiki.thebiogrid.org/doku.php/biogridrest
api_key_file = os.path.dirname(os.path.realpath(__file__)) + '/' + \
               'biogrid_api_key'
api_key_env_name = 'BIOGRID_API_KEY'

# Try to read the API key from a file
try:
    with open(api_key_file, 'rt') as fh:
        api_key = fh.read().strip()
except IOError:
    logger.error('BioGRID API key could not be found, trying environment '
                 'variable $%s.' % api_key_env_name)
    logger.error(api_key_file)
    # Try the environment variable
    if api_key_env_name in os.environ:
        api_key = os.environ.get(api_key_env_name)
    else:
        logger.error('No BioGRID API key found in environment variable '
                     '%s.' % api_key_env_name)
        api_key = None


def get_interactors(gene_name):
    res_dict = _send_request([gene_name], include_interactors=True)
    interaction_list = []
    for result in res_dict.values():
        if result['OFFICIAL_SYMBOL_A'] == gene_name and \
           result['OFFICIAL_SYMBOL_B'] == gene_name:
            interaction_list.append(gene_name)
        elif result['OFFICIAL_SYMBOL_A'] == gene_name:
            interaction_list.append(result['OFFICIAL_SYMBOL_B'])
        elif result['OFFICIAL_SYMBOL_B'] == gene_name:
            interaction_list.append(result['OFFICIAL_SYMBOL_A'])
        else:
            assert False, "Interaction doesn't contain target gene!"
    interaction_counter = Counter(interaction_list)
    interaction_counter = sorted(interaction_counter.items(),
                                 key=lambda x: x[1], reverse=True)
    return interaction_counter


def get_statements(gene_list):
    res_dict = _send_request(gene_list, include_interactors=True)
    statements = []
    for int_id, interaction in res_dict.items():
        agent_a_name = interaction['OFFICIAL_SYMBOL_A']
        agent_b_name = interaction['OFFICIAL_SYMBOL_B']
        agent_a = Agent(agent_a_name, db_refs={'HGNC': agent_a_name})
        agent_b = Agent(agent_b_name, db_refs={'HGNC': agent_b_name})
        ev = Evidence(source_api='biogrid',
                      source_id=int_id,
                      pmid=interaction['PUBMED_ID'],
                      text=None,
                      annotations=interaction)
        stmt = Complex([agent_a, agent_b], evidence=ev)
        statements.append(stmt)
    return statements


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
        # The json module produces strings, not bytes, so the file should be
        # opened in text mode
        with open(save_json_name, 'wt') as fh:
            json.dump(res_dict, fh, indent=1)
    publications = _extract_publications(res_dict, gene_names)
    return publications


@python_2_unicode_compatible
class Publication(object):
    def __init__(self, interaction, interaction_id):
        self.pmid = "PMID" + str(interaction['PUBMED_ID'])
        self.modification = interaction['MODIFICATION']
        self.experimental_system = interaction['EXPERIMENTAL_SYSTEM']
        self.experimental_system_type = interaction['EXPERIMENTAL_SYSTEM_TYPE']
        self.throughput = interaction['THROUGHPUT']
        self.interaction_id = interaction_id

    def __str__(self):
        return "Publication(%s)" % self.pmid

    def __repr__(self):
        return str(self)


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


def _send_request(gene_names, include_interactors=False):
    if api_key is None:
        logger.error('BioGRID cannot be used without API key')
        return None
    params = {'searchNames': 'true',
              'geneList': '|'.join(gene_names),
              'taxId': '9606',
              'format': 'json',
              'includeInteractors': include_interactors,
              'accesskey': api_key}
    res = requests.get(biogrid_url, params)
    res.raise_for_status()
    # The json module handles the conversion from bytes to unicode internally
    res_dict = res.json()
    return res_dict
