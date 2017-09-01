from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import ndex.client
from .processor import NdexCxProcessor

def process_cx_file(file_name):
    with open(file_name, 'rt') as fh:
        json_list = json.load(fh)
        return process_cx(json_list)

def process_ndex_network(network_id, username=None, password=None):
    nd = ndex.client.Ndex(username=username, password=password)
    res = nd.get_network_as_cx_stream(network_id)
    if res.status_code != 200:
        logger.error('Problem downloading network: status code %s' %
                     res.status_code)
        logger.error('Response: %s' % res.text)
        return None
    json_list = res.json()
    return process_cx(json_list)

def process_cx(json_list):
    ncp = NdexCxProcessor(json_list)
    ncp.get_statements()
    return ncp
