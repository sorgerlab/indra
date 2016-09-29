from __future__ import print_function, unicode_literals
import xml.etree.ElementTree as et
from indra.util import UnicodeXMLTreeBuilder as UTB
from indra.util import unicode_strs

try:
    from io import BytesIO as SIO
except ImportError:
    from StringIO import StringIO as SIO

def test_unicode_tree_builder():
    xml = u'<html><bar>asdf</bar></html>'.encode('utf-8')
    xml_io = SIO(xml)
    tree = et.parse(xml_io, parser=UTB())
    bar = tree.find('.//bar')
    assert unicode_strs(bar)

