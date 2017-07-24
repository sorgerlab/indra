from __future__ import unicode_literals
from builtins import dict, str
import logging
import requests

logger = logging.getLogger("omnipath")

_loaded = True

op_url = 'http://omnipathdb.org/'

def get_ptms(gene_list):
    params = {'format': 'json'}
    gene_str = ','.join(gene_list)
    ptm_url = '%s/ptms/%s' % (op_url, gene_str)
    res = requests.get(ptm_url, params=params)
    if not res.status_code == 200 or not res.text:
        return None
    return res.json


