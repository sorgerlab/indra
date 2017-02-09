from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import xml.etree.ElementTree as ET
from indra.util import UnicodeXMLTreeBuilder as UTB
from indra.sparser.processor import SparserProcessor

logger = logging.getLogger('sparser')

def process_xml(xml_str):
    try:
        tree = ET.XML(xml_str, parser=UTB())
    except ET.ParseError as e:
        logger.error('Could not parse XML string')
        logger.error(e)
        return None
    sp = _process_elementtree(tree)
    return sp


def _process_elementtree(tree):
    sp = SparserProcessor(tree)
    sp.get_modifications()
    sp.get_activations()
    return sp
