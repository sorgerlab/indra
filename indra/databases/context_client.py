import json
import requests

ndexbio_context = 'http://general.bigmech.ndexbio.org:8081/context/'

def get_protein_expression(gene_names, cell_types):
    req_type = 'expression/cell_line'
    if isinstance(gene_names, basestring):
        gene_names = [gene_names]
    if isinstance(cell_types, basestring):
        cell_types = [cell_types]
    params = {g: cell_types for g in gene_names}
    res = send_request(req_type, params)
    return res

def get_mutations(gene_names, cell_types):
    req_type = 'mutation/cell_line'
    if isinstance(gene_names, basestring):
        gene_names = [gene_names]
    if isinstance(cell_types, basestring):
        cell_types = [cell_types]
    params = {g: cell_types for g in gene_names}
    res = send_request(req_type, params)
    return res

def send_request(req_type, params=None):
    if params is None:
        params = {}
    res = requests.post(ndexbio_context + req_type, json=params)
    if res.status_code != 200:
        print 'Request to NDEx service returned with status %d' % res.status_code
        return None
    res_json = res.json()
    return res_json
