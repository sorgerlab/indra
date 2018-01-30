from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import json
import xml.etree.ElementTree as ET
from indra.util import UnicodeXMLTreeBuilder as UTB, kappy_json_to_graph
from indra.util import unicode_strs
from io import BytesIO
from os import path


def test_unicode_tree_builder():
    xml = u'<html><bar>asdf</bar></html>'.encode('utf-8')
    xml_io = BytesIO(xml)
    tree = ET.parse(xml_io, parser=UTB())
    bar = tree.find('.//bar')
    assert unicode_strs(bar)


def test_kappy_influence_json_to_graph():
    with open(path.join(path.dirname(path.abspath(__file__)), 'kappy_influence.json'), 'r') as f:
        imap = json.load(f)
    graph = kappy_json_to_graph(imap)
    assert graph is not None, 'No graph produced.'
    n_nodes = len(graph.nodes)
    n_edges = len(graph.edges)
    assert n_nodes == 13, \
        'Wrong number (%d vs. %d) of nodes on the graph.' % (n_nodes, 13)
    assert n_edges == 6, \
        "Wrong number (%d vs. %d) of edges on graph." % (n_edges, 4)

