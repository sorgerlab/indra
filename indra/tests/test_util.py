from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from io import BytesIO
import xml.etree.ElementTree as ET
from indra.util import unicode_strs
from indra.util import UnicodeXMLTreeBuilder as UTB


def test_unicode_tree_builder():
    xml = u'<html><bar>asdf</bar></html>'.encode('utf-8')
    xml_io = BytesIO(xml)
    tree = ET.parse(xml_io, parser=UTB())
    bar = tree.find('.//bar')
    assert unicode_strs(bar)
