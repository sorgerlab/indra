from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers import PybelAssembler
import indra.tools.assemble_corpus as ac
import pybel

def assemble_pybel(stmts, out_file_prefix):
    """Return a PyBEL Assembler"""
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    pba = PybelAssembler(stmts)
    pba.make_model()
    with open(out_file_prefix, 'wt') as f:
        pybel.to_json_file(pba.model, f)
