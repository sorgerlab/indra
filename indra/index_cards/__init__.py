from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
from indra.index_cards.processor import IndexCardProcessor

def process_json_file(file_name, source_api):
    with open(file_name, 'rt') as fh:
        json_dict = json.load(fh)
        return process_json(json_dict, source_api)

def process_json_str(json_str, source_api):
    json_dict = json.loads(json_str)
    return process_json(json_dict, source_api)

def process_json(json_dict, source_api):
    if isinstance(json_dict, dict):
        json_dict = [json_dict]
    icp = IndexCardProcessor(json_dict, source_api)
    icp.get_modifications()
    icp.get_complexes()
    icp.get_binds()
    icp.get_translocates()
    return icp
