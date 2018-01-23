import os
from collections import Counter
import pygraphviz
import networkx as nx
import numpy as np
from indra.explanation.paths_graph import *

source = 'A'
target = 'D'
target_polarity = 0

graph1_s = nx.DiGraph()
graph1_s.add_nodes_from(['A', 'B', 'C', 'D'])
graph1_s.add_edges_from([('A', 'B', {'sign': 0}),
                         ('B', 'D', {'sign': 0}),
                         ('A', 'C', {'sign': 0}),
                         ('C', 'D', {'sign': 0})])

graph1_uns = nx.DiGraph()
graph1_uns.add_nodes_from(['A', 'B', 'C', 'D'])
graph1_uns.add_edges_from([('A', 'B'), ('B', 'D'), ('A', 'C'), ('C', 'D')])


g_samp = nx.DiGraph() # Graph for testing sampling uniformly vs. non-uniformly
g_samp.add_edges_from([
    ('source', 'A1'), ('source', 'A2'),
    ('A1', 'B1'),
    ('A2', 'B2'), ('A2', 'B3'), ('A2', 'B4'), ('A2', 'B5'),
    ('B1', 'target'),
    ('B2', 'target'), ('B3', 'target'), ('B4', 'target'), ('B5', 'target')])


g_split_nodes = nx.DiGraph()
g_split_nodes.add_edges_from(
    [(0, 1), (0, 2), (0, 3), (0, 5), (1, 0), (1, 2), (1, 3), (1, 4), (2, 0),
     (2, 3), (2, 4), (3, 0), (3, 1), (3, 2), (3, 4), (3, 5), (4, 1), (4, 2),
     (4, 3), (5, 1), (5, 3), (5, 4)])


def test_get_reachable_sets_unsigned():
    f_level, b_level = get_reachable_sets(graph1_uns, source, target,
                                          signed=False)
    assert f_level == {0: {'A'}, 1: {'B', 'C'}, 2: {'D'}}
    assert b_level == {0: {'D'}, 1: {'B', 'C'}, 2: {'A'}}


def test_get_reachable_sets_signed():
    f_level, b_level = get_reachable_sets(graph1_s, source, target, signed=True)
    assert f_level == {0: {('A', 0)}, 1: {('B', 0), ('C', 0)}, 2: {('D', 0)}}
    assert b_level == {0: {('D', 0)}, 1: {('B', 0), ('C', 0)}, 2: {('A', 0)}}


def test_unreachability_unsigned():
    graph = nx.DiGraph()
    graph.add_nodes_from(['A', 'B', 'C', 'D'])
    graph.add_edges_from([('A', 'B'), ('D', 'B'), ('C', 'A'), ('C', 'D')])
    (f_level, b_level) = get_reachable_sets(graph, source, target,
                                            max_depth=5, signed=False)
    assert f_level == {}
    assert b_level == {}


def test_unreachability_signed():
    # First make the unreachability due to the direction of the edges
    graph = nx.DiGraph()
    graph.add_nodes_from(['A', 'B', 'C', 'D'])
    graph.add_edges_from([('A', 'B', {'sign': 0}),
                          ('D', 'B', {'sign': 0}),
                          ('C', 'A', {'sign': 0}),
                          ('C', 'D', {'sign': 0})])
    (f_level, b_level) = get_reachable_sets(graph, source, target,
                                            max_depth=5, signed=True)
    assert f_level == {}
    assert b_level == {}
    # This time, make the unreachability due to the sign
    graph = nx.DiGraph()
    graph.add_nodes_from(['A', 'B', 'C', 'D'])
    graph.add_edges_from([('A', 'B', {'sign': 0}),
                          ('D', 'B', {'sign': 0}),
                          ('C', 'A', {'sign': 0}),
                          ('C', 'D', {'sign': 0})])
    (f_level, b_level) = get_reachable_sets(graph, source, target, max_depth=5,
                                            signed=True)
    assert f_level == {}
    assert b_level == {}


def test_from_graph_unsigned():
    # Path length 1
    f_level, b_level = get_reachable_sets(graph1_s, source, target, max_depth=3,
                                          signed=False)
    pg = PathsGraph.from_graph(graph1_uns, source, target, 1, f_level, b_level,
                               signed=False)
    assert len(pg.graph) == 0
    # Path length 2
    pg = PathsGraph.from_graph(graph1_uns, source, target, 2, f_level, b_level, 
                               signed=False)
    paths = list(nx.shortest_simple_paths(pg.graph, (0, 'A'), (2, 'D')))
    assert len(paths) == 2
    assert [(0, 'A'), (1, 'C'), (2, 'D')] in paths
    assert [(0, 'A'), (1, 'B'), (2, 'D')] in paths
    # Path length 3
    pg = PathsGraph.from_graph(graph1_uns, source, target, 3, f_level, b_level, 
                               signed=False)
    assert len(pg.graph) == 0


def test_from_graph_unsigned_no_levels():
    length = 2
    pg = PathsGraph.from_graph(graph1_uns, source, target, length)
    assert isinstance(pg, PathsGraph)
    paths = list(nx.shortest_simple_paths(pg.graph, (0, 'A'), (2, 'D')))
    assert len(paths) == 2
    assert [(0, 'A'), (1, 'C'), (2, 'D')] in paths
    assert [(0, 'A'), (1, 'B'), (2, 'D')] in paths


def test_from_graph_signed():
    # Path length 1
    f_level, b_level = get_reachable_sets(graph1_s, source, target, signed=True,
                                          max_depth=3)
    pg = PathsGraph.from_graph(graph1_s, source, target, 1, f_level, b_level,
                               signed=True, target_polarity=0)
    assert not pg.graph
    # Path length 2
    pg = PathsGraph.from_graph(graph1_s, source, target, 2, f_level, b_level,
                               signed=True, target_polarity=0)
    paths = list(nx.shortest_simple_paths(pg.graph, (0, ('A', 0)),
                                                          (2, ('D', 0))))
    assert len(paths) == 2
    assert [(0, ('A', 0)), (1, ('C', 0)), (2, ('D', 0))] in paths
    assert [(0, ('A', 0)), (1, ('B', 0)), (2, ('D', 0))] in paths
    # Path length 3
    pg = PathsGraph.from_graph(graph1_s, source, target, 3, f_level, b_level,
                               signed=True, target_polarity=0)
    assert not pg.graph


def test_pg_check_unreachable_unsigned():
    graph = nx.DiGraph()
    graph.add_nodes_from(['A', 'B', 'C', 'D'])
    graph.add_edges_from([('A', 'B'), ('D', 'B'), ('C', 'A'), ('C', 'D')])
    (f_level, b_level) = get_reachable_sets(graph, source, target, max_depth=5, 
                                            signed=False)
    assert f_level == {}
    assert b_level == {}
    pg = PathsGraph.from_graph(graph, source, target, 2, f_level, b_level,
                               signed=False)
    assert not pg.graph
    # A graph where there is a path, but not of the given length (3)
    (f_level, b_level) = get_reachable_sets(graph1_s, source, target,
                                            max_depth=5, signed=False)
    pg = PathsGraph.from_graph(graph, source, target, 3, f_level, b_level,
                               signed=False)
    assert not pg.graph


def test_multidigraph_signed():
    graph = nx.MultiDiGraph()
    graph.add_edges_from([('A', 'B', {'sign': 0}), ('A', 'B', {'sign': 1})])
    f_level, b_level = get_reachable_sets(graph, 'A', 'B', max_depth=3,
                                          signed=True)
    assert f_level[0] == {('A', 0)}
    assert f_level[1] == {('B', 0), ('B', 1)}
    assert b_level[0] == {('B', 0)}
    assert b_level[1] == {('A', 0), ('A', 1)}


def test_sample_paths():
    g_uns = nx.DiGraph()
    g_uns.add_edges_from((('A', 'B'), ('A', 'C'), ('C', 'D'), ('B', 'D'),
                          ('D', 'B'), ('D', 'C'), ('B', 'E'), ('C', 'E')))
    source, target, length = ('A', 'E', 4)
    pg = PathsGraph.from_graph(g_uns, source, target, length)
    sample_paths = pg.sample_paths(100)
    assert set(sample_paths) == set(
        [('A', 'B', 'D', 'B', 'E'),
         ('A', 'B', 'D', 'C', 'E'),
         ('A', 'C', 'D', 'B', 'E'),
         ('A', 'C', 'D', 'C', 'E')])


def test_enumerate_paths():
    g_uns = nx.DiGraph()
    g_uns.add_edges_from((('A', 'B'), ('A', 'C'), ('C', 'D'), ('B', 'D'),
                          ('D', 'B'), ('D', 'C'), ('B', 'E'), ('C', 'E')))
    source, target, length = ('A', 'E', 4)
    pg = PathsGraph.from_graph(g_uns, source, target, length)
    enum_paths = pg.enumerate_paths()
    assert set(enum_paths) == set(
        [('A', 'B', 'D', 'B', 'E'),
         ('A', 'B', 'D', 'C', 'E'),
         ('A', 'C', 'D', 'B', 'E'),
         ('A', 'C', 'D', 'C', 'E')])


def test_count_paths():
    g_uns = nx.DiGraph()
    g_uns.add_edges_from((('A', 'B'), ('A', 'C'), ('C', 'D'), ('B', 'D'),
                          ('D', 'B'), ('D', 'C'), ('B', 'E'), ('C', 'E')))
    source, target, length = ('A', 'E', 4)
    pg = PathsGraph.from_graph(g_uns, source, target, length)
    num_paths = pg.count_paths()
    assert num_paths == 4


def test_non_uniform_sampling():
    pg = PathsGraph.from_graph(g_samp, 'source', 'target', 3)
    # There are five different paths, but sampling uniformly based on local
    # edge weights should result in ~50% of paths going through B1
    os.environ['TEST_FLAG'] = 'TRUE'
    np.random.seed(1) # Seed the random number generator
    num_samples = 1000
    paths = pg.sample_paths(num_samples)
    num_b1_paths = len([p for p in paths if 'B1' in p])
    num_other_paths = len([p for p in paths if 'B1' not in p])
    assert num_b1_paths == 526
    assert num_other_paths == 474


def test_uniform_sampling():
    pg = PathsGraph.from_graph(g_samp, 'source', 'target', 3)
    # There are five different paths; sampling uniformly across the whole
    # path distribution should result in 20% of paths going through each of
    # paths going through B1-B5.
    pg.set_uniform_path_distribution()
    os.environ['TEST_FLAG'] = 'TRUE'
    np.random.seed(1) # Seed the random number generator
    num_samples = 5000
    path_count = pg.count_paths()
    assert path_count == 5
    paths = pg.sample_paths(num_samples, weighted=True)
    b_ctr = Counter([p[2] for p in paths])
    print(b_ctr)
    assert b_ctr == {'B1': 1021, 'B2': 991, 'B3': 964, 'B4': 1022, 'B5': 1002}


def test_raw_graph_sampling():
    """Smoke test for sampling on an underlying graph."""
    g = nx.DiGraph()
    g.add_edges_from([('S', 'A'), ('S', 'T'), ('A', 'T'), ('A', 'S')])
    paths = sample_raw_graph(g, 'S', 'T', num_samples=10000,
                             eliminate_cycles=False)
    path_ctr = Counter([len(p) for p in paths])
    print(path_ctr)

    paths = sample_raw_graph(g, 'S', 'T', num_samples=10000,
                             eliminate_cycles=True)
    path_ctr = Counter([len(p) for p in paths])
    assert len(path_ctr) == 2
    print(path_ctr)

    paths = sample_raw_graph(g, 'S', 'T', num_samples=10000,
                             eliminate_cycles=False, max_depth=4)
    path_ctr = Counter([len(p) for p in paths])
    assert len(path_ctr) == 4
    print(path_ctr)


def test_combine_paths_graphs():
    g = nx.DiGraph()
    g.add_edges_from([('S', 'A'), ('S', 'T'), ('A', 'T'), ('A', 'S')])
    max_depth = 4
    pg_dict = {}
    for length in range(1, max_depth+1):
        paths_graph = PathsGraph.from_graph(g, 'S', 'T', length)
        pg_dict[length] = paths_graph
    cpg = combine_paths_graphs(pg_dict)
    paths = cpg.sample_paths(1000)
    path_ctr = Counter(paths)


def test_combine_cfpgs():
    g = nx.DiGraph()
    g.add_edges_from([('S', 'A'), ('S', 'T'), ('A', 'T'), ('A', 'S')])
    max_depth = 4
    pg_dict = {}
    for length in range(1, max_depth+1):
        cfpg = CFPG.from_graph(g, 'S', 'T', length)
        pg_dict[length] = cfpg
    cpg = combine_cfpgs(pg_dict)
    paths = cpg.sample_paths(1000)
    path_ctr = Counter(paths)


def test_combine_cfpgs_link_split_nodes():
    g = g_split_nodes
    # Make combined paths graph
    pg_dict = {}
    max_depth = 5
    src, tgt = (0, 5)
    for length in range(max_depth+1):
        cfpg = CFPG.from_graph(g, src, tgt, length)
        pg_dict[length] = cfpg
    cpg = combine_cfpgs(pg_dict)
    # The edges that should be added by the edge addition procedure:
    linked_edges = [
       ((2, 4, frozenset({(0, 0), (1, 1), (2, 4)})),
        (3, 3, frozenset({(0, 0), (1, 1), (1, 2), (2, 2), (2, 4), (3, 3)}))),
       ((2, 4, frozenset({(0, 0), (1, 2), (2, 4)})),
        (3, 3, frozenset({(0, 0), (1, 1), (1, 2), (2, 2), (2, 4), (3, 3)}))),
     ]
    for e in linked_edges:
        assert e in cpg.graph.edges()


def test_correct_edge_multiplicities():
    g = g_split_nodes
    # Make combined paths graph
    pg_dict = {}
    max_depth = 5
    src, tgt = (0, 5)
    for length in range(max_depth+1):
        cfpg = CFPG.from_graph(g, src, tgt, length)
        pg_dict[length] = cfpg
    cpg = combine_cfpgs(pg_dict)
    cpg.correct_edge_multiplicities()
    # The edges that should be adjusted, each with 0.5
    adjust_edges = [
       ((1, 2, frozenset({(1, 2), (0, 0)})),
        (2, 4, frozenset({(1, 2), (0, 0), (1, 1), (2, 4)}))),
       ((1, 2, frozenset({(1, 2), (0, 0)})),
        (2, 4, frozenset({(1, 2), (0, 0), (2, 4)}))),
       ((1, 1, frozenset({(0, 0), (1, 1)})),
        (2, 4, frozenset({(1, 2), (0, 0), (1, 1), (2, 4)}))),
       ((1, 1, frozenset({(0, 0), (1, 1)})),
        (2, 4, frozenset({(0, 0), (1, 1), (2, 4)})))]
    for u, v, data in cpg.graph.edges_iter(data=True):
        assert 'weight' in data
        if (u, v) in adjust_edges:
            assert data['weight'] == 0.5
        else:
            assert data['weight'] == 1


def test_path_tree():
    g = nx.DiGraph()
    g.add_edges_from((('A', 'B'), ('A', 'C'), ('A', 'E'),
                          ('B', 'D'), ('B', 'E'),
                          ('C', 'D'), ('C', 'E'),
                          ('D', 'B'), ('D', 'C'), ('D', 'E'),
                          ))
    source, target, length = ('A', 'E', 4)
    paths = list(nx.all_simple_paths(g, source, target))
    pt = PathTree(paths)
    pt_ref_edges = set([
            (tuple(), ('A',)),
            (('A',), ('A', 'E')),
            (('A',), ('A', 'B')),
            (('A',), ('A', 'C')),
            (('A', 'B'), ('A', 'B', 'E')),
            (('A', 'B'), ('A', 'B', 'D')),
            (('A', 'C'), ('A', 'C', 'E')),
            (('A', 'C'), ('A', 'C', 'D')),
            (('A', 'B', 'D'), ('A', 'B', 'D', 'E')),
            (('A', 'C', 'D'), ('A', 'C', 'D', 'E')),
            (('A', 'B', 'D'), ('A', 'B', 'D', 'C')),
            (('A', 'C', 'D'), ('A', 'C', 'D', 'B')),
            (('A', 'B', 'D', 'C'), ('A', 'B', 'D', 'C', 'E')),
            (('A', 'C', 'D', 'B'), ('A', 'C', 'D', 'B', 'E')),
        ])
    assert set(pt.graph.edges()) == pt_ref_edges

    # Sample from the path tree
    num_samples = 1000
    samp_paths = pt.sample(num_samples=num_samples)
    assert len(samp_paths) == num_samples
    assert set(samp_paths) == set([tuple(p) for p in paths])


def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')


if __name__ == '__main__':
    edge_prob = 0.5
    num_nodes = 10
    src = 0
    tgt = num_nodes - 1
    max_depth = num_nodes - 1
    num_samples = 30000

    #while True:
    #    g = nx.random_graphs.erdos_renyi_graph(num_nodes, edge_prob,
    #                                               directed=True)
    g = nx.DiGraph()
    g.add_edges_from([
 (0, 1),
 (0, 2),
 (0, 3),
 (0, 5),
 (0, 8),
 (1, 3),
 (1, 4),
 (1, 5),
 (1, 7),
 (1, 8),
 (1, 9),
 (2, 4),
 (2, 5),
 (2, 7),
 (2, 8),
 (2, 9),
 (3, 1),
 (3, 5),
 (3, 6),
 (3, 7),
 (3, 8),
 (3, 9),
 (4, 2),
 (4, 3),
 (4, 5),
 (4, 8),
 (5, 0),
 (5, 2),
 (5, 3),
 (5, 4),
 (5, 9),
 (6, 3),
 (6, 4),
 (6, 5),
 (6, 8),
 (6, 9),
 (7, 0),
 (7, 2),
 (7, 5),
 (7, 8),
 (8, 1),
 (8, 7),
 (8, 9),
 (9, 0),
 (9, 1),
 (9, 2),
 (9, 4),
 (9, 6)])

    # Make combined paths graph
    pg_list = []
    for length in range(max_depth+1):
        cfpg = CFPG.from_graph(g, src, tgt, length)
        pg_list.append(cfpg)
    cpg = CombinedCFPG(pg_list)

    # Build path tree
    all_paths = list(nx.all_simple_paths(g, src, tgt))
    pt = PathTree(all_paths)

    # Draw graphs
    #draw(g, 'problem_graph.pdf')
    #draw(pt.graph, 'problem_graph_tree.pdf')
    #draw(cpg.graph, 'problem_combined_cfpg.pdf')

    print("Sampling tree")
    tree_paths = pt.sample(num_samples)
    tree_path_ctr = Counter(tree_paths)

    print("Sampling CFPG")
    cfpg_paths = cpg.sample_paths(num_samples)
    cfpg_path_ctr = Counter(cfpg_paths)
    assert set(tree_path_ctr.keys()) == set(cfpg_path_ctr.keys())

    # Build contingency table for statistical comparison
    tree_path_ctr = sorted([(k, v) for k, v in tree_path_ctr.items()],
                            key=lambda x: x[0])
    cfpg_path_ctr = sorted([(k, v) for k, v in cfpg_path_ctr.items()],
                           key=lambda x: x[0])
    # Check to make sure we have paths; if not, both should be empty
    #if cfpg_path_ctr == [] and tree_path_ctr == []:
    #    continue

    cont_tbl = [[t[1] for t in cfpg_path_ctr],
                [t[1] for t in tree_path_ctr]]
    import scipy.stats
    chi2_result = scipy.stats.chi2_contingency(cont_tbl)
    pval = chi2_result[1]
    print("Tree")
    print(tree_path_ctr)
    print(Counter(tree_paths))
    print("CFPG")
    print(cfpg_path_ctr)
    print(Counter(cfpg_paths))
    print(chi2_result)
    #assert pval > 0.01



"""


g.add_edges_from(
[(0, 1),
(0, 2),
(0, 3),
(0, 5),
(1, 0),
(1, 2),
(1, 3),
(1, 4),
(2, 0),
(2, 3),
(2, 4),
(3, 0),
(3, 1),
(3, 2),
(3, 4),
(3, 5),
(4, 1),
(4, 2),
(4, 3),
(5, 1),
(5, 3),
(5, 4)])

g = nx.DiGraph()
g.add_edges_from([(0, 1),
 (0, 3),
 (1, 4),
 (1, 5),
 (2, 0),
 (2, 1),
 (2, 5),
 (3, 0),
 (3, 1),
 (3, 4),
 (4, 0),
 (4, 1),
 (4, 2),
 (4, 3),
 (5, 2),
 (5, 3),
 (5, 4)])


"""
