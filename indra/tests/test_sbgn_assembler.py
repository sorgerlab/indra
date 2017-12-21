from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import requests
from lxml import etree
from indra.statements import *
from indra.assemblers.sbgn_assembler import SBGNAssembler

ns = {'sbgn': 'http://sbgn.org/libsbgn/pd/0.1'}
schema_url = 'https://raw.githubusercontent.com/sbgn/libsbgn/master/resources/SBGN.xsd'

def _get_parser():
    res = requests.get(schema_url)
    xsd_str = res.text
    schema_root = etree.XML(xsd_str)
    schema = etree.XMLSchema(schema_root)
    # FIXME: the validation should be reinstated when it works
    #parser = etree.XMLParser(schema=schema)
    parser = etree.XMLParser()
    return parser

sbgn_parser = _get_parser()

def _parse_sbgn(sbgn_xml):
    print(sbgn_xml)
    et = etree.fromstring(sbgn_xml, sbgn_parser)
    return et

def _test_numelements(et, nglyphs, narcs):
    glyphs = et.findall('sbgn:map/sbgn:glyph', namespaces=ns)
    assert(len(glyphs) == nglyphs)
    arcs = et.findall('sbgn:map/sbgn:arc', namespaces=ns)
    assert(len(arcs) == narcs)

def test_modification():
    st = Phosphorylation(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    et = _parse_sbgn(sbgn_xml)
    _test_numelements(et, 4, 3)

def test_remove_modification():
    st = Deacetylation(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    et = _parse_sbgn(sbgn_xml)
    _test_numelements(et, 4, 3)

def test_activation():
    st = Activation(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    et = _parse_sbgn(sbgn_xml)
    _test_numelements(et, 4, 3)

def test_inhibition():
    st = Inhibition(Agent('BRAF'), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    et = _parse_sbgn(sbgn_xml)
    _test_numelements(et, 4, 3)

def test_increaseamount():
    st = IncreaseAmount(Agent(''), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    et = _parse_sbgn(sbgn_xml)
    _test_numelements(et, 3, 2)

def test_bound_condition():
    bc = BoundCondition(Agent('RAF1'), True)
    st = Phosphorylation(Agent('BRAF', bound_conditions=[bc]), Agent('MAP2K1'))
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    et = _parse_sbgn(sbgn_xml)
    _test_numelements(et, 4, 3)

def test_complex():
    egfr1 = Agent('EGFR',
                  bound_conditions=[BoundCondition(Agent('EGF'), True)])
    egfr2 = Agent('EGFR',
                  bound_conditions=[BoundCondition(Agent('EGF'), True)])
    st = Complex([egfr1, egfr2])
    sa = SBGNAssembler([st])
    sbgn_xml = sa.make_model()
    et = _parse_sbgn(sbgn_xml)
    _test_numelements(et, 4, 3)

def test_activeform():
    erkact = Agent('MAPK1', activity=ActivityCondition('kinase', True))
    erkelk = Phosphorylation(erkact, Agent('ELK1'))
    erkp = Agent('MAPK1', mods=[ModCondition('phosphorylation')])
    st = ActiveForm(erkp, 'kinase', True)
    sa = SBGNAssembler([st, erkelk])
    sbgn_xml = sa.make_model()
    et = _parse_sbgn(sbgn_xml)
    _test_numelements(et, 4, 3)
