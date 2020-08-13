import numpy as np
import networkx as nx

from indra.explanation.pathfinding.pathfinding import bfs_search, \
    shortest_simple_paths, bfs_search_multiple_nodes
from indra.explanation.model_checker.model_checker import \
    signed_edges_to_signed_nodes

INT_PLUS = 0
INT_MINUS = 1


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
    edge_by_hash = {
        'HASH1': [
            ('A1', 'B1'),
            ('B1', 'C1'),
            ('B3', 'C1'),
        ],
        'HASH2': [
            ('B2', 'C1'),
            ('B1', 'C1'),
            ('A2', 'B1'),
            ('A3', 'B2'),
            ('A4', 'B2'),
        ]
    }

    return edges, signed_edges, edge_beliefs, list(all_ns), edge_by_hash


def _setup_unsigned_graph():
    edges, signed_edges, edge_beliefs, all_ns, edge_by_hash = _digraph_setup()
    dg = nx.DiGraph()
    dg.add_edges_from(edges)

    # Add belief
    for e in dg.edges:
        dg.edges[e]['belief'] = edge_beliefs[e]
        dg.edges[e]['weight'] = -np.log(edge_beliefs[e], dtype=np.longfloat)
    
    # Add edge_by_hash
    dg.graph['edge_by_hash'] = edge_by_hash

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
    edges, signed_edges, edge_beliefs, all_ns, edge_by_hash = _digraph_setup()
    seg = nx.MultiDiGraph()

    seg.add_edges_from(signed_edges)
    # ATTN!! seg.edges yields u, v, index while seg.edges() yields u, v
    for u, v, sign in seg.edges:
        seg.edges[(u, v, sign)]['sign'] = sign
        seg.edges[(u, v, sign)]['belief'] = edge_beliefs[(u, v)]

    for node, data in seg.nodes(data=True):
        data['ns'] = node[0]
        data['id'] = node[1]

    sng = signed_edges_to_signed_nodes(graph=seg, prune_nodes=True,
                                       copy_edge_data=True)
    for u, v in sng.edges:
        sng.edges[(u, v)]['weight'] = -np.log(sng.edges[(u, v)]['belief'])
    
    seg.graph['edge_by_hash'] = edge_by_hash
    sng.graph['edge_by_hash'] = edge_by_hash

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
    expected_paths = {('D1', 'C1', 'B1'), ('D1', 'C1', 'B2'),
                      ('D1', 'C1', 'B3')}
    paths = [p for p in bfs_search(dg, 'D1', depth_limit=5,
                                   reverse=True, terminal_ns=['b'],
                                   node_filter=all_ns)]
    assert len(paths) == len(expected_paths), len(paths)
    assert set(paths) == expected_paths, 'sets of paths not equal'

    # Terminate on 'a'
    expected_paths = {('D1', 'C1', 'B1', 'A1'), ('D1', 'C1', 'B1', 'A2'),
                      ('D1', 'C1', 'B2', 'A3'), ('D1', 'C1', 'B2', 'A4')}
    paths = [p for p in bfs_search(dg, 'D1', depth_limit=5,
                                   reverse=True, terminal_ns=['a'],
                                   node_filter=all_ns)]
    assert len(paths) == len(expected_paths), len(paths)
    assert set(paths) == expected_paths, 'sets of paths not equal'

    # Test memory limit; a very low number should yield one path
    error = False
    gen = bfs_search(dg, 'D1', depth_limit=5, reverse=True, max_memory=16)
    _ = next(gen)
    try:
        _ = next(gen)
    except (RuntimeError, StopIteration):
        error = True
    assert error


def test_signed_bfs():
    seg, sng, all_ns = _setup_signed_graph()
    # D1 being upregulated: 13 paths
    paths = [p for p in bfs_search(
        g=sng, source_node=('D1', INT_PLUS), reverse=True, depth_limit=5,
        node_filter=all_ns, sign=INT_PLUS)
    ]
    assert len(paths) == 13, len(paths)


def test_shortest_simple_paths_mod_unsigned():
    dg, all_ns = _setup_unsigned_graph()
    dg.add_edge('B1', 'A3', belief=0.7, weight=-np.log(0.7))  # Add long path
    # between
    # B1 and C1
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

    # Weighted searches
    paths = [p for p in shortest_simple_paths(dg, source, target,
                                              weight='weight')]
    assert tuple(paths[0]) == ('B1', 'C1', 'D1')
    assert tuple(paths[1]) == ('B1', 'A3', 'B2', 'C1', 'D1')


def test_shortest_simple_paths_mod_signed():
    seg, sng, all_ns = _setup_signed_graph()

    # Add a very long path
    sng.add_edge(('B3', INT_PLUS), ('Z1', INT_PLUS),
                 belief=0.6, weight=-np.log(0.6))

    source = ('B2', INT_PLUS)
    target = ('D1', INT_PLUS)  # D1 upregulated

    # Check that beliafs and weights have been set
    assert sng.edges[(source, ('C1', INT_PLUS))].get('weight', False)
    assert sng.edges[(source, ('C1', INT_PLUS))].get('belief', False)

    # Unweighted search
    expected_paths = {(source, ('C1', INT_PLUS), target),
                      (source, ('C1', INT_MINUS), target)}
    paths = [tuple(p) for p in shortest_simple_paths(sng, source, target)]
    assert len(paths) == 2
    assert set(paths) == expected_paths, 'sets of paths not equal'

    # Weighted search
    source = ('B3', INT_PLUS)
    expected_paths = {
        (source, ('C1', 1), target),
        (source, ('Z1', 0), ('A1', 0), ('B1', 0), ('C1', 0), target),
        (source, ('Z1', 0), ('A1', 1), ('B1', 1), ('C1', 1), target)
    }
    paths = tuple([tuple(p) for p in shortest_simple_paths(sng, source, target,
                                                           weight='weight')])
    assert len(paths) == 3
    assert set(paths) == expected_paths, 'sets of paths not equal'


def test_bfs_multiple_nodes():
    dg, all_ns = _setup_unsigned_graph()

    # Run regular bfs search for nodes separately
    paths = [p for p in bfs_search(dg, 'C1', depth_limit=2, reverse=True)]
    assert len(paths) == 7, len(paths)
    paths = [p for p in bfs_search(dg, 'D1', depth_limit=2, reverse=True)]
    assert len(paths) == 4, len(paths)
    # number of paths with several nodes bfs should be a sum of separate runs
    paths = [p for p in bfs_search_multiple_nodes(
        dg, ['C1', 'D1'], depth_limit=2, reverse=True)]
    assert len(paths) == 11, len(paths)
    # Path limit should limit the number of paths (both across all paths and
    # within one node search)
    paths = [p for p in bfs_search_multiple_nodes(
        dg, ['C1', 'D1'], depth_limit=2, reverse=True, path_limit=9)]
    assert len(paths) == 9, len(paths)
    paths = [p for p in bfs_search_multiple_nodes(
        dg, ['C1', 'D1'], depth_limit=2, reverse=True, path_limit=5)]
    print(len(paths))
    assert len(paths) == 5, len(paths)


def test_shortest_simple_paths_strict_mesh_id_filtering():
    G = _setup_unsigned_graph()[0]
    G.add_edge('A2', 'B3', belief=0.7, weight=-np.log(0.7))
    G.add_edge('B3', 'B1', belief=0.7, weight=-np.log(0.7))
    G.graph['edge_by_hash']['HASH1'] += [('A2', 'B3'), ('B3', 'B1')]
    
    def count_paths(source, target, hashes):
        try:
            paths = [p for p in shortest_simple_paths(G, source, target, 
                                                      hashes=hashes, 
                                                      weight='weight',
                                                      strict_mesh_id_filtering=True)]
            return len(paths)
        except nx.NetworkXNoPath:
            return 0

    assert count_paths('Z1', 'C1', []) == 1
    assert count_paths('A1', 'C1', []) == 1
    assert count_paths('A2', 'C1', []) == 3
    assert count_paths('A2', 'D1', []) == 3

    assert count_paths('Z1', 'C1', ['HASH1']) == 0
    assert count_paths('A1', 'C1', ['HASH1']) == 1
    assert count_paths('A2', 'C1', ['HASH1']) == 2
    assert count_paths('A2', 'D1', ['HASH1']) == 0
    
    assert count_paths('Z1', 'C1', ['HASH2']) == 0
    assert count_paths('A1', 'C1', ['HASH2']) == 0
    assert count_paths('A2', 'C1', ['HASH2']) == 1
    assert count_paths('A2', 'D1', ['HASH2']) == 0

    assert count_paths('Z1', 'C1', ['HASH1', 'HASH2']) == 0
    assert count_paths('A1', 'C1', ['HASH1', 'HASH2']) == 1
    assert count_paths('A2', 'C1', ['HASH1', 'HASH2']) == 3
    assert count_paths('A2', 'D1', ['HASH1', 'HASH2']) == 0


def test_shortest_simple_paths_weighed_by_mesh_ids():
    G = _setup_unsigned_graph()[0]
    G.add_edge('A3', 'B1', belief=0.7, weight=-np.log(0.7))
    G.graph['edge_by_hash']['HASH1'].append(('A3', 'B1'))
    source = 'A3'
    target = 'C1'
    paths = list(shortest_simple_paths(G, source, target, hashes=['HASH1']))
    assert paths == [['A3', 'B1', 'C1'], ['A3', 'B2', 'C1']]
    paths = list(shortest_simple_paths(G, source, target, hashes=['HASH2']))
    assert paths == [['A3', 'B2', 'C1'], ['A3', 'B1', 'C1']]
    paths = list(shortest_simple_paths(G, source, target, hashes=['HASH1', 'HASH2']))
    assert paths == [['A3', 'B1', 'C1'], ['A3', 'B2', 'C1']]


def test_bfs_strict_mesh_id_filtering():
    dg = _setup_unsigned_graph()[0]

    paths = [p for p in bfs_search(dg, 'C1', depth_limit=6, reverse=True, 
                                   strict_mesh_id_filtering=True, hashes=[])]
                                   
    expected = {('C1', 'B3'), 
                ('C1', 'B2'), 
                ('C1', 'B2', 'A3'), 
                ('C1', 'B2', 'A4'), 
                ('C1', 'B1'), 
                ('C1', 'B1', 'A2'),
                ('C1', 'B1', 'A1'),
                ('C1', 'B1', 'A1', 'Z1')}
    assert len(paths) == 8
    assert set(paths) == expected
    
    paths = [p for p in bfs_search(dg, 'C1', depth_limit=6, reverse=True, 
                                   strict_mesh_id_filtering=True, hashes=['HASH1'])]
    expected = {('C1', 'B3'),
                ('C1', 'B1'),
                ('C1', 'B1', 'A1')}
    assert len(paths) == 3
    assert set(paths) == expected

    paths = [p for p in bfs_search(dg, 'C1', depth_limit=6, reverse=True,
                                   strict_mesh_id_filtering=True, hashes=['HASH2'])]
    expected = {('C1', 'B2'),
                ('C1', 'B2', 'A3'),
                ('C1', 'B2', 'A4'),
                ('C1', 'B1'),
                ('C1', 'B1', 'A2')}
    assert len(paths) == 5
    assert set(paths) == expected

    paths = [p for p in bfs_search(dg, 'C1', depth_limit=6, reverse=True,
                                   strict_mesh_id_filtering=True,
                                   hashes=['HASH1', 'HASH2'])]
    expected = {('C1', 'B3'),
                ('C1', 'B2'),
                ('C1', 'B2', 'A3'),
                ('C1', 'B2', 'A4'),
                ('C1', 'B1'),
                ('C1', 'B1', 'A2'),
                ('C1', 'B1', 'A1')}
    assert len(paths) == 7
    assert set(paths) == expected
