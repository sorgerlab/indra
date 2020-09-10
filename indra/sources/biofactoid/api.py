import requests
from .processor import BioFactoidProcessor

biofactoid_url = 'https://biofactoid.org/api/document'
biofactoid_unstable_url = 'https://unstable.factoid.baderlab.org/api/document'


def process_from_web(url=None):
    url = url if url else biofactoid_url
    res = requests.get(url)
    res.raise_for_status()
    return process_json(res.json())


def process_json(biofactoid_json):
    bp = BioFactoidProcessor(biofactoid_json)
    bp.extract_statements()
    return bp