from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
from .processor import NdexCxProcessor

def process_cx_file(file_name):
    with open(file_name, 'rt') as fh:
        json_dict = json.load(fh)
        return process_cx(json_dict)

def process_cx_str(json_str):
    json_dict = json.loads(json_str)
    return process_cx(json_dict)

def process_cx(json_dict):
    if isinstance(json_dict, dict):
        json_dict = [json_dict]
    ncp = NdexCxProcessor(json_dict)
    ncp.get_modifications()
    ncp.get_complexes()
    ncp.get_binds()
    return ncp
