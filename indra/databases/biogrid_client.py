import os
import json
import logging
import requests

biogrid_url = 'http://webservice.thebiogrid.org/interactions/'

logger = logging.getLogger('biopax')

# THIS FILE IS NOT UNDER VERSION CONTROL
# For more information see http://wiki.thebiogrid.org/doku.php/biogridrest
api_key_file = os.path.dirname(os.path.realpath(__file__)) + '/' + \
               'biogrid_api_key'
# Read the API key
try:
    with open(api_key_file, 'rt') as fh:
        api_key = fh.read().strip()
except IOError:
    logging.error('BioGRID API key could not be found.')
    logging.error(api_key_file)
    api_key = None

def process_query(statement, save_json_name='biogrid_output.json'):
    agent_names = [agent.name for agent in statement.agent_list() \
                    if agent is not None]
    if len(agent_names) < 2:
        return None
    res_dict = _send_request(agent_names)
    with open(save_json_name, 'wt') as fh:
        json.dump(res_dict, fh, indent=1)
    return process_json(res_dict, agent_names)

def get_biogrid_interactors(dict):
    b_list = [dict['OFFICIAL_SYMBOL_A'], dict['OFFICIAL_SYMBOL_B']]
    return b_list

def get_subset(dict, agent_list):
    new_dict = {}
    for id in dict.keys():
        b_list = get_biogrid_interactors(dict[id])
        if set(b_list) == set(agent_list):
            new_dict[id] = dict[id]
    return new_dict

def process_json(dict, agent_list):
    agent_subset = get_subset(dict, agent_list)
    publications = []
    for paper in agent_subset.keys():
        publications.append(Publication(agent_subset[paper]))
    return publications

def _send_request(agent_list):
    params = {'searchNames': 'true',
              'geneList': '|'.join(agent_list),
              'taxId': '9606',
              'format': 'jsonExtended',
              'accesskey': api_key}
    res = requests.get(biogrid_url, params)
    res.raise_for_status()
    res_dict = res.json()
    return res_dict

class Publication(object):
    def __init__(self, interaction):
        self.pmid = "PMID" + str(interaction['PUBMED_ID'])
        self.modification = interaction['MODIFICATION']
        self.experimental_system = interaction[
            'EXPERIMENTAL_SYSTEM']
        self.experimental_system_type = interaction[
            'EXPERIMENTAL_SYSTEM_TYPE']
        self.throughput = interaction['THROUGHPUT']


    def __str__(self):
        return ("PMID:" + self.pmid)


    def __repr__(self):
        return ("<pmid: %s at 0x%x>" %
                (self.pmid, id(self)))
