__all__ = ['process_xml_file']

from xml.etree import ElementTree as ET
from .processor import SemRepXmlProcessor


def process_xml_file(fname, use_gilda_grounding=False, predicate_mappings=None):
    tree = ET.parse(fname)
    sp = SemRepXmlProcessor(tree, use_gilda_grounding=use_gilda_grounding,
                            predicate_mappings=predicate_mappings)
    return sp
