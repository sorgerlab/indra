from __future__ import print_function
from indra.util import unicode_parser
import xml.etree.ElementTree as et

xml = u'<html><bar>asdf</bar></html>'.encode('utf-8')
try:
    from io import BytesIO
    xml_io = BytesIO(xml)
except ImportError:
    from StringIO import StringIO
    xml_io = StringIO(xml)

tree = et.parse(xml_io, parser=unicode_parser)
bar = tree.find('.//bar')
print(type(bar.text))

