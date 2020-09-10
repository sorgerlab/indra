import requests
from .processor import BioFactoidProcessor

biofactoid_url = 'https://biofactoid.org/api/document'


def process_from_web():
    res = requests.get(biofactoid_url)
    res.raise_for_status()
    return process_json(res.json())


def process_json(biofactoid_json):
    bp = BioFactoidProcessor(biofactoid_json)
    bp.extract_statements()
    return bp