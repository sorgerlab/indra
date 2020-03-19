import requests
from indra.config import get_config
from .processor import HypothesisProcessor


base_url = 'https://api.hypothes.is/api/'
api_key = get_config('HYPOTHESIS_API_KEY')
headers = {'Authorization': 'Bearer %s' % api_key,
           'Accept': 'application/vnd.hypothesis.v1+json',
           'content-type': 'application/json'}


def send_request(endpoint, **params):
    if api_key is None:
        return ValueError('No API key set in HYPOTHESIS_API_KEY')
    res = requests.get(base_url + endpoint, headers=headers,
                       params=params)
    res.raise_for_status()
    return res.json()


def process_annotations(group):
    res = send_request('search', group=group)
    annotations = res.get('rows', [])
    hp = HypothesisProcessor(annotations)
    hp.get_statements()
    return hp
