from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import xml.etree.ElementTree as ET
from indra.statements import *
from indra.assemblers.sbgn_assembler import SBGNAssembler

ns = {'sbgn': 'http://sbgn.org/libsbgn/pd/0.1'}

def _test_numelements(sbgn_xml, nglyphs, narcs):
    et = ET.fromstring(sbgn_xml)
    glyphs = et.findall('sbgn:map/sbgn:glyph', namespaces=ns)
    assert(len(glyphs) == nglyphs)
    arcs = et.findall('sbgn:map/sbgn:arc', namespaces=ns)
    assert(len(arcs) == narcs)

def test_modification():
    st = Phosphorylation(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    _test_numelements(sbgn_xml, 4, 3)

def test_remove_modification():
    st = Deacetylation(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    _test_numelements(sbgn_xml, 4, 3)

def test_activation():
    st = Activation(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    _test_numelements(sbgn_xml, 4, 3)

def test_inhibition():
    st = Inhibition(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    _test_numelements(sbgn_xml, 4, 3)

def test_increaseamoutn():
    st = IncreaseAmount(Agent(''), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    _test_numelements(sbgn_xml, 4, 3)

