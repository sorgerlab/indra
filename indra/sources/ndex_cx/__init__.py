from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
from .processor import NdexCxProcessor

def process_cx_file(file_name):
    with open(file_name, 'rt') as fh:
        json_list = json.load(fh)
        return process_cx(json_list)

def process_cx_str(json_str):
    json_list = json.loads(json_str)
    return process_cx(json_list)

def process_cx(json_list):
    ncp = NdexCxProcessor(json_list)
    ncp.get_statements()
    return ncp
