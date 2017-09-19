from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers import PybelAssembler
import indra.tools.assemble_corpus as ac
import pybel
import requests

def assemble_pybel(stmts, out_file_prefix):
    """Return a PyBEL Assembler"""
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)

    pba = PybelAssembler(stmts, name='INDRA/REACH Korkut Model',
                         description='Automatically assembled model of '
                                     'cancer signaling.',
                         version='0.0.10')
    pba.make_model()
    pybel.to_bel_path(pba.model, out_file_prefix + '.bel')
    with open(out_file_prefix, 'wt') as f:
        pybel.to_json_file(pba.model, f)
    url =  'https://pybel.scai.fraunhofer.de/api/receive'
    headers = {'content-type': 'application/json'}
    requests.post(url, json=pybel.to_json(pba.model), headers=headers)
