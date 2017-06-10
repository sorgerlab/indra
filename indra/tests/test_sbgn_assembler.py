from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.assemblers.sbgn_assembler import SBGNAssembler

def test_modification():
    st = Phosphorylation(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    print(sbgn_xml)
