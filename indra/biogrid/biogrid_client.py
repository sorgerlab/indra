import requests
import json

for line in open('authentication.txt').readlines():
    access_key = line.strip()
    
base_url = 'http://webservice.thebiogrid.org/interactions/'


def get_json(agent_list):

    params = {'searchNames': 'true',
              'geneList': '|'.join(agent_list),
              'taxId': '9606',
              'format': 'jsonExtended',
              'accesskey': access_key}
              
    r = requests.get(base_url, params)
    r.raise_for_status()
    dict = r.json()
    return dict


def save_json(dict, filename):
    with open(filename, 'w') as f:
        json.dump(dict, f)
