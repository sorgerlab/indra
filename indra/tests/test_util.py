from __future__ import print_function, unicode_literals
import xml.etree.ElementTree as ET
from indra.util import UnicodeXMLTreeBuilder as UTB
from indra.util import unicode_strs
from io import BytesIO

def test_unicode_tree_builder():
    xml = u'<html><bar>asdf</bar></html>'.encode('utf-8')
    xml_io = BytesIO(xml)
    tree = ET.parse(xml_io, parser=UTB())
    bar = tree.find('.//bar')
    assert unicode_strs(bar)

