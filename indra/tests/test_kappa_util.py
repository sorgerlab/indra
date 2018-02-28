from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import json
from os import path
from indra.util.kappa_util import im_json_to_graph, cm_json_to_graph


def test_kappy_influence_json_to_graph():
    with open(path.join(path.dirname(path.abspath(__file__)),
                        'kappy_influence.json'), 'r') as f:
        imap = json.load(f)
    graph = im_json_to_graph(imap)
    assert graph is not None, 'No graph produced.'
    n_nodes = len(graph.nodes())
    n_edges = len(graph.edges())
    assert n_nodes == 13, \
        'Wrong number (%d vs. %d) of nodes on the graph.' % (n_nodes, 13)
    assert n_edges == 6, \
        "Wrong number (%d vs. %d) of edges on graph." % (n_edges, 4)


def test_kappy_contact_json_to_graph():
    with open(path.join(path.dirname(path.abspath(__file__)),
                        'kappy_contact.json'), 'r') as f:
        cmap = json.load(f)
    graph = cm_json_to_graph(cmap)
    assert graph is not None, 'No graph produced.'
    n_nodes = len(graph.nodes())
    n_edges = len(graph.edges())
    n_subgraphs = len(graph.subgraphs())
    assert n_nodes == 6, \
        'Wrong number (%d vs. %d) of nodes on the graph.' % (n_nodes, 6)
    assert n_edges == 3, \
        "Wrong number (%d vs. %d) of edges on graph." % (n_edges, 3)
    assert n_subgraphs == 4, \
        "Wrong number (%d vs. %d) of subgraphs on graph." % (n_subgraphs, 4)
