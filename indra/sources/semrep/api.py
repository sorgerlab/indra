__all__ = ['process_xml_file']

from xml.etree import ElementTree as ET
from .processor import SemRepXmlProcessor


def process_xml_file(fname):
    tree = ET.parse(fname)
    sp = SemRepXmlProcessor(tree)
    return sp