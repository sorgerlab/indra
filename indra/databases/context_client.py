import json
import requests

ndexbio_services = 'http://general.bigmech.ndexbio.org:8081/context/'

def get_protein_expression(gene_names, cell_types):
    req_type = 'expression/cell_line'
    params = {g: cell_types for g in gene_names}
    res = send_request(req_type, params)
    return res

def get_mutations(gene_names, cell_types):
    req_type = 'mutation/cell_line'
    params = {g: cell_types for g in gene_names}
    res = send_request(req_type, params)
    return res

def send_request(req_type, params):
    res = requests.post(ndexbio_services + req_type, json=params)
    res_json = res.json()
    return res_json
