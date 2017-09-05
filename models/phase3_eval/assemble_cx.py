from indra.util import _require_python3
from os.path import join as pjoin
from indra.assemblers import CxAssembler
import indra.tools.assemble_corpus as ac

def assemble_cx(stmts, out_file_prefix, network_type):
    """Return a CX assembler."""
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    if network_type == 'direct':
        stmts = ac.filter_direct(stmts)

    out_file = '%s_%s.cx' % (out_file_prefix, network_type)

    ca = CxAssembler()
    ca.add_statements(stmts)
    model = ca.make_model()
    ca.save_model(out_file)
    return ca
