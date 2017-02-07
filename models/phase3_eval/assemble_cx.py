from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join as pjoin
from indra.assemblers import CxAssembler
import indra.tools.assemble_corpus as ac

def assemble_cx(stmts, out_file):
    """Return a CX assembler."""
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)

    ca = CxAssembler()
    ca.add_statements(stmts)
    model = ca.make_model()
    ca.save_model(out_file)
    return ca
