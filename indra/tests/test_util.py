from io import BytesIO
import xml.etree.ElementTree as ET

from indra.statements import Statement

from indra.util import unicode_strs
from indra.util import UnicodeXMLTreeBuilder as UTB
from indra.util.statement_presentation import _get_relation_keyed_stmts


def test_unicode_tree_builder():
    xml = u'<html><bar>asdf</bar></html>'.encode('utf-8')
    xml_io = BytesIO(xml)
    tree = ET.parse(xml_io, parser=UTB())
    bar = tree.find('.//bar')
    assert unicode_strs(bar)


def test_conversion_keying():
    stmt_json = {"type": "Conversion",
                 "subj": {"name": "inflammatory response", "db_refs": {}},
                 "obj_from": [{"name": "KNG1",
                               "db_refs": {"HGNC": "6383", "UP": "P01042"}}],
                 "obj_to": [{"name": "Kallidin",
                             "db_refs": {"SCHEM": "Kallidin"}}],
                 "id": "d2361669-dfe5-45e0-914a-c96a82ad25fb"}
    stmt_list = [Statement._from_json(stmt_json)]
    stmt_list[0].agent_list()
    list(_get_relation_keyed_stmts(stmt_list))
    return

