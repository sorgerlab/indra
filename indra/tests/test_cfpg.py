import os
import pickle
from collections import Counter
from os.path import dirname, join
import numpy as np
import networkx as nx
from indra.explanation import paths_graph as pg

random_graph_pkl = join(dirname(__file__), 'random_graphs.pkl')

g_uns = nx.DiGraph()
g_uns.add_edges_from((('A', 'B'), ('A', 'C'), ('C', 'D'), ('B', 'D'),
                      ('D', 'B'), ('D', 'C'), ('B', 'E'), ('C', 'E')))
source = 'A'
target = 'E'
length = 4


def test_from_graph_no_levels():
    cfpg = pg.CFPG.from_graph(g_uns, source, target, length)
    assert isinstance(cfpg, pg.CFPG)
    paths = cfpg.enumerate_paths()
    assert len(paths) == 2
    assert ('A', 'B', 'D', 'C', 'E') in paths
    assert ('A', 'C', 'D', 'B', 'E') in paths
    assert len(cfpg.graph) == 8
    # The D node should be split into two nodes
    d_nodes = [n for n in cfpg.graph.nodes() if n[1] == 'D']
    assert len(d_nodes) == 2


def test_from_graph_with_levels_bad_depth():
    """Raise an exception if the requested path length is greater than the
    depth of the provided reach sets."""
    (f_reach, b_reach) = pg.get_reachable_sets(g_uns, source, target,
                                                        max_depth=2)
    cfpg = pg.CFPG.from_graph(g_uns, source, target, length,
                              fwd_reachset=f_reach, back_reachset=b_reach)
    assert not cfpg.graph


def test_from_pg():
    (f_reach, b_reach) = pg.get_reachable_sets(g_uns, source, target,
                                                        max_depth=length)
    pg_0 = pg.PathsGraph.from_graph(g_uns, source, target, length, f_reach,
                                      b_reach)
    cfpg = pg.CFPG.from_pg(pg_0)
    paths = cfpg.enumerate_paths()
    assert len(paths) == 2
    assert ('A', 'B', 'D', 'C', 'E') in paths
    assert ('A', 'C', 'D', 'B', 'E') in paths
    assert len(cfpg.graph) == 8
    # The D node should be split into two nodes
    d_nodes = [n for n in cfpg.graph.nodes() if n[1] == 'D']
    assert len(d_nodes) == 2


def test_sample_paths():
    cfpg = pg.CFPG.from_graph(g_uns, source, target, length)
    sample_paths = cfpg.sample_paths(100)
    assert set(sample_paths) == set(
        [('A', 'B', 'D', 'C', 'E'),
         ('A', 'C', 'D', 'B', 'E')])


def test_enumerate_paths():
    cfpg = pg.CFPG.from_graph(g_uns, source, target, length)
    enum_paths = cfpg.enumerate_paths()
    assert set(enum_paths) == set(
        [('A', 'B', 'D', 'C', 'E'),
         ('A', 'C', 'D', 'B', 'E')])


def test_count_paths():
    cfpg = pg.CFPG.from_graph(g_uns, source, target, length)
    num_paths = cfpg.count_paths()
    assert num_paths == 2


def test_on_random_graphs():
    """For each of 25 random graphs, check that the number of cycle free paths
    for a given depth and source/target pair matches the results from
    networkx all_simple_paths. Graphs range from rough"""
    # We use 25 randomly generated graphs for testing the algorithm
    with open(random_graph_pkl, 'rb') as f:
        rg_dict = pickle.load(f)

    min_depth = 5
    max_depth = 10
    for i in range(1):
        G_i, source, target = rg_dict[i]
        print("graph# %d, %d nodes, %d edges" % (i, len(G_i), len(G_i.edges())))
        (f_reach, b_reach)  = pg.get_reachable_sets(G_i, source, target,
                                        max_depth=max_depth, signed=False)
        # Try different path lengths
        for length in range(min_depth, max_depth+1):
            print("Checking paths of length %d" % length)
            # For validation, we compute explicitly the set of paths in the
            # original graph of a fixed length
            P = list(nx.all_simple_paths(G_i, source, target, length+1))
            # Filter to paths of this length
            P_correct = [tuple(p) for p in P if len(p) == length+1]
            # Generate the raw paths graph
            G_cf = pg.CFPG.from_graph(G_i, source, target, length, f_reach,
                                      b_reach)
            # Check the path count
            path_count = G_cf.count_paths()
            assert len(P_correct) == path_count
            # Enumerate paths using node tuples
            P_cf_pruned = G_cf.enumerate_paths(names_only=False)
            # Next we extract the paths by projecting down to second
            # component (node names)
            P_cf_pruned_names = G_cf.enumerate_paths(names_only=True)
            print("# of paths: %d" % len(P_cf_pruned_names))

            # We verify the three required properties.
            # Recall:
            # CF1: Every source-to-target path in G_cf is cycle free.
            # CF2: Every cycle free path in the original graph appears as a
            #      source-to-target path in G_cf.
            # CF3: There is a 1-1 correspondence between the paths in G_cf and
            # the paths in the original graph. This means there is no
            # redundancy in the representation. For every path in the original
            # graph there is a unique path in G_cf that corresponds to it.

            # We first verify CF1.
            for p in P_cf_pruned_names:
                if len(p) != len(list(set(p))):
                    print("cycle!")
                    print(p)
                    assert False
            # Next we verify CF2. We will in fact check if the set of paths in
            # P_cf_pruned_names is exactly the set of paths in the original
            # graph.
            if set(P_correct) != set(P_cf_pruned_names):
                print("Paths do not match reference set from networkx")
                print("graph, length", (i, length))
                assert False
            # Finally we verify CF3
            if len(P_cf_pruned) != len(list(set(P_cf_pruned_names))):
                print("redundant representation!")
                print("graph, length", (i, length))
                assert False


def test_example_graph1():
    sif_file = join(dirname(__file__), 'korkut_im.sif')
    g = pg.load_signed_sif(sif_file)
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'
    target_polarity = 0
    enum_paths = pg.enumerate_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=False)
    assert len(enum_paths) == len(set(enum_paths))
    enum_paths = set(enum_paths)
    sampled_paths = set(pg.sample_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8,
                             num_samples=10000, cycle_free=False))
    count_paths = pg.count_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=False)
    assert len(sampled_paths) == len(enum_paths)
    assert len(sampled_paths) == count_paths
    assert sampled_paths == enum_paths

    enum_cf_paths = pg.enumerate_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=True)
    assert len(enum_cf_paths) == len(set(enum_cf_paths))
    enum_cf_paths = set(enum_cf_paths)
    sampled_cf_paths = set(pg.sample_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8,
                             num_samples=10000, cycle_free=True))
    count_cf_paths = pg.count_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=True)
    assert len(sampled_cf_paths) == len(enum_cf_paths)
    assert len(sampled_cf_paths) == count_cf_paths
    assert sampled_cf_paths == enum_cf_paths

    # 101 cycle free paths
    assert count_cf_paths == 101
    # All of the cycle-free paths should be contained in the set of paths
    # with cycles
    assert sampled_paths.intersection(sampled_cf_paths) == sampled_cf_paths
    # Check that all paths in the set difference contain cycles
    for path in sampled_paths.difference(sampled_cf_paths):
        assert len(set(path)) < len(path)


def test_uniform_sampling_example_graph1():
    sif_file = join(dirname(__file__), 'korkut_im.sif')
    g = pg.load_signed_sif(sif_file)
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'
    target_polarity = 0
    length = 8
    cfpg = pg.CFPG.from_graph(g, source, target, length, signed=True,
                              target_polarity=0)
    os.environ['TEST_FLAG'] = 'TRUE'
    np.random.seed(1)
    # Count paths
    # Now, re-weight for uniformity and re-sample
    num_samples = cfpg.count_paths() * 1000
    cfpg.set_uniform_path_distribution()
    sampled_paths_uni = cfpg.sample_paths(num_samples=num_samples)
    ctr_uni = Counter(sampled_paths_uni)
    for path, count in ctr_uni.items():
        assert count > 900 and count < 1100


def test_enumerate_example_graph2():
    sif_file = join(dirname(__file__), 'korkut_stmts.sif')
    g = pg.load_signed_sif(sif_file)
    source = 'BLK'
    target = 'EIF4EBP1'
    enum_paths = pg.enumerate_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=False)
    assert len(enum_paths) == len(set(enum_paths))
    enum_paths = set(enum_paths)
    count_paths = pg.count_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=False)

    enum_cf_paths = pg.enumerate_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=True)
    assert len(enum_cf_paths) == len(set(enum_cf_paths))
    enum_cf_paths = set(enum_cf_paths)
    count_cf_paths = pg.count_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=True)
    # Check that the counts match the enumeratio
    assert len(enum_paths) == count_paths
    assert len(enum_cf_paths) == count_cf_paths
    # All of the cycle-free paths should be contained in the set of paths
    # with cycles
    assert enum_paths.intersection(enum_cf_paths) == enum_cf_paths
    # Check that all paths in the set difference contain cycles
    for path in enum_paths.difference(enum_cf_paths):
        assert len(set(path)) < len(path)


def test_sampling_example_graph2():
    sif_file = join(dirname(__file__), 'korkut_stmts.sif')
    g = pg.load_signed_sif(sif_file)
    source = 'BLK'
    target = 'EIF4EBP1'
    enum_cf_paths = set(pg.enumerate_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=True))
    sampled_cf_paths = set(pg.sample_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8,
                             num_samples=10000, cycle_free=True))
    count_cf_paths = pg.count_paths(g, source, target, signed=True,
                             target_polarity=0, max_depth=8, cycle_free=True)
    # The sampled paths should be a subset of the enumerated paths
    assert sampled_cf_paths.intersection(enum_cf_paths) == sampled_cf_paths


def test_weighted_sampling():
    g = nx.DiGraph()
    g.add_edges_from([
        ('A', 'B', {'weight': 3}),
        ('A', 'C', {'weight': 1}),
        ('C', 'D'),
        ('B', 'D'),
        ('D', 'B'),
        ('D', 'C'),
        ('B', 'E'),
        ('C', 'E')])
    source, target, length = ('A', 'E', 4)
    cfpg = pg.CFPG.from_graph(g, source, target, length)
    os.environ['TEST_FLAG'] = 'TRUE'
    np.random.seed(1)
    samp_paths = cfpg.sample_paths(1000)
    ctr = Counter(samp_paths)
    assert ctr[('A', 'B', 'D', 'C', 'E')] == 767
    assert ctr[('A', 'C', 'D', 'B', 'E')] == 233


def test_combine_cfpgs():
    g = nx.DiGraph()
    g.add_edges_from([('S', 'A'), ('S', 'T'), ('A', 'T'), ('A', 'S')])
    max_depth = 4
    pg_list = []
    for length in range(1, max_depth+1):
        cfpg = pg.CFPG.from_graph(g, 'S', 'T', length)
        pg_list.append(cfpg)
    cpg = pg.CombinedCFPG(pg_list)
    paths = cpg.sample_paths(1000)
    path_ctr = Counter(paths)


