import numpy as np
import networkx as nx

from indra.explanation.pathfinding.pathfinding import bfs_search, \
    shortest_simple_paths
from indra.explanation.pathfinding.util import signed_edges_to_signed_nodes

INT_PLUS, INT_MINUS = 0, 1


def _digraph_setup():
    # Ensures alphabetical order in reverse traversal
    edge_beliefs = {('Z1', 'A1'): 1 - 0.2,
                    ('A1', 'B1'): 1 - 0.2,
                    ('A2', 'B1'): 1 - 0.3,
                    ('A3', 'B2'): 1 - 0.5,
                    ('A4', 'B2'): 1 - 0.6,
                    ('B1', 'C1'): 1 - 0.2,
                    ('B2', 'C1'): 1 - 0.3,
                    ('B3', 'C1'): 1 - 0.4,
                    ('C1', 'D1'): 1 - 0.2}
    edges = [('Z1', 'A1'),
             ('A1', 'B1'),
             ('A2', 'B1'),
             ('A3', 'B2'),
             ('A4', 'B2'),
             ('B1', 'C1'),
             ('B2', 'C1'),
             ('B3', 'C1'),
             ('C1', 'D1')]
    signed_edges = [
        ('Z1', 'A1', INT_PLUS),   # 1
        ('Z1', 'A1', INT_MINUS),  # 2
        ('A1', 'B1', INT_PLUS),   # 3
        ('A2', 'B1', INT_MINUS),  # 4
        ('B1', 'C1', INT_PLUS),   # 5
        ('A3', 'B2', INT_PLUS),   # 6
        ('A4', 'B2', INT_MINUS),  # 7
        ('B2', 'C1', INT_PLUS),   # 8
        ('B2', 'C1', INT_MINUS),  # 9
        ('B3', 'C1', INT_MINUS),  # 10
        ('C1', 'D1', INT_PLUS),   # 11
        ('C1', 'D1', INT_MINUS),  # 12
    ]
    all_ns = set()
    for e in edges:
        all_ns.add(e[0][0].lower())
        all_ns.add(e[1][0].lower())

    return edges, signed_edges, edge_beliefs, list(all_ns)


def _setup_unsigned_graph():
    edges, signed_edges, edge_beliefs, all_ns = _digraph_setup()
    dg = nx.DiGraph()
    dg.add_edges_from(edges)

    # Add belief
    for e in dg.edges:
        dg.edges[e]['belief'] = edge_beliefs[e]
        dg.edges[e]['weight'] = -np.log(dg.edges[e]['belief'],
                                        dtype=np.longfloat)

    # Add namespaces
    nodes1, nodes2 = list(zip(*edges))
    nodes = set(nodes1).union(nodes2)
    for node in nodes:
        ns = node[0]
        _id = node[1]
        dg.nodes[node]['ns'] = ns
        dg.nodes[node]['id'] = _id
    return dg, all_ns


def _setup_signed_graph():
    edges, signed_edges, edge_beliefs, all_ns = _digraph_setup()
    seg = nx.MultiDiGraph()

    seg.add_edges_from(signed_edges)
    # ATTN!! seg.edges yields u, v, index while seg.edges() yields u, v
    for u, v, sign in seg.edges:
        seg.edges[(u, v, sign)]['sign'] = sign
        seg.edges[(u, v, sign)]['belief'] = edge_beliefs[(u, v)]

    sng = signed_edges_to_signed_nodes(graph=seg, prune_nodes=True,
                                       copy_edge_data=False)
    return seg, sng, all_ns


def test_bfs():
    dg, all_ns = _setup_unsigned_graph()

    # Test basic part of algorithm
    paths = [p for p in bfs_search(dg, 'C1', depth_limit=1, reverse=True)]
    assert len(paths) == 3, len(paths)
    paths = [p for p in bfs_search(dg, 'C1', depth_limit=2, reverse=True)]
    assert len(paths) == 7, len(paths)
    paths = [p for p in bfs_search(dg, 'C1', depth_limit=2, reverse=True,
                                   path_limit=4)]
    assert len(paths) == 4, len(paths)

    # Test ns allowance list
    ans = ['c', 'b']
    assert len([p for p in
                bfs_search(dg, 'C1', depth_limit=2, reverse=True,
                           node_filter=ans)]) == 3
    assert all(len(p) < 3 for p in
               bfs_search(dg, 'C1', depth_limit=2, reverse=True,
                          node_filter=ans))

    # Test longer paths
    assert len([p for p in bfs_search(dg, 'D1', depth_limit=5,
                                      reverse=True)]) == 9

    # Test node blacklist
    assert len([p for p in bfs_search(dg, 'D1', depth_limit=5, reverse=True,
                                      node_blacklist={'Z1'})]) == 8

    # Test max per node option
    # Should get 4 paths with max_per_node=1
    expected_paths = {('D1', 'C1'), ('D1', 'C1', 'B1'),
                      ('D1', 'C1', 'B1', 'A1'),
                      ('D1', 'C1', 'B1', 'A1', 'Z1')}
    paths = [p for p in bfs_search(dg, 'D1', depth_limit=5,
                                   reverse=True, max_per_node=1,
                                   node_filter=all_ns)]
    assert len(paths) == 4, len(paths)
    assert set(paths) == expected_paths, 'sets of paths not equal'

    # Test terminal NS
    # Terminate on 'b'
    expected_paths = {('D1', 'C1'), ('D1', 'C1', 'B1'), ('D1', 'C1', 'B2'),
                      ('D1', 'C1', 'B3')}
    paths = [p for p in bfs_search(dg, 'D1', depth_limit=5,
                                   reverse=True, terminal_ns=['b'],
                                   node_filter=all_ns)]
    assert len(paths) == 4, len(paths)
    assert set(paths) == expected_paths, 'sets of paths not equal'
    # Terminate on 'a'
    expected_paths = {('D1', 'C1'), ('D1', 'C1', 'B1'), ('D1', 'C1', 'B2'),
                      ('D1', 'C1', 'B3'), ('D1', 'C1', 'B1', 'A1'),
                      ('D1', 'C1', 'B1', 'A2'), ('D1', 'C1', 'B2', 'A3'),
                      ('D1', 'C1', 'B2', 'A4')}
    paths = [p for p in bfs_search(dg, 'D1', depth_limit=5,
                                   reverse=True, terminal_ns=['a'],
                                   node_filter=all_ns)]
    assert len(paths) == len(expected_paths), len(paths)
    assert set(paths) == expected_paths, 'sets of paths not equal'


def test_signed_bfs():
    dg, _ = _setup_unsigned_graph()
    seg, sng, all_ns = _setup_signed_graph()
    # D1 being upregulated: 12 paths
    paths = [p for p in bfs_search(
        g=sng, source_node=('D1', INT_PLUS), g_nodes=dg.nodes,
        g_edges=seg.edges, reverse=True, depth_limit=5, node_filter=all_ns,
        sign=INT_PLUS)
    ]
    assert len(paths) == 13, len(paths)


def test_shortest_simple_paths_mod_unsigned():
    dg, all_ns = _setup_unsigned_graph()
    dg.add_edge('B1', 'A3', belief=0.7)  # Create long path between B1 and C1
    source, target = 'B1', 'D1'

    # Unweighted searches
    paths = [p for p in shortest_simple_paths(dg, source, target)]
    assert len(paths) == 2
    assert tuple(paths[0]) == ('B1', 'C1', 'D1')
    assert tuple(paths[1]) == ('B1', 'A3', 'B2', 'C1', 'D1')
    assert len([p for p in shortest_simple_paths(
        dg, source, target, ignore_nodes={'A3'})]) == 1
    # Test nx.NoPathFound
    try:
        len([p for p in shortest_simple_paths(
            dg, source, target, ignore_nodes={'C1'})]) == 0
    except Exception as exc:
        assert isinstance(exc, nx.NetworkXNoPath)

    # Weigthed searches
    paths = [p for p in shortest_simple_paths(dg, source, target,
                                              weight='weight')]
    assert tuple(paths[0]) == ('B1', 'C1', 'D1')
    assert tuple(paths[1]) == ('B1', 'A3', 'B2', 'C1', 'D1')


def test_shortest_simple_paths_mod_signed():
    seg, sng, all_ns = _setup_signed_graph()
    source = ('B2', INT_PLUS)
    target = ('D1', INT_PLUS)  # D1 upregulated
    expected_paths = {(source, ('C1', INT_PLUS), target),
                      (source, ('C1', INT_MINUS), target)}
    paths = [tuple(p) for p in shortest_simple_paths(sng, source, target)]
    assert len(paths) == 2
    assert set(paths) == expected_paths, 'sets of paths not equal'
