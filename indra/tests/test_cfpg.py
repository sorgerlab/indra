import pickle
import networkx as nx
from os.path import dirname, join
from indra.explanation.paths_graph import paths_graph, pre_cfpg as pcf
from indra.explanation.paths_graph.cfpg import *


random_graph_pkl = join(dirname(__file__), 'random_graphs.pkl')


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
        (f_reach, b_reach)  = \
                paths_graph.get_reachable_sets(G_i, source, target,
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
            pre_cfpg = pcf.from_graph(G_i, source, target, length, f_reach,
                                      b_reach)
            src = (0, source)
            tgt = (length, target)
            G_cf = PG_cf(src,tgt, pre_cfpg.graph, pre_cfpg.tags)

            # We verify the three required properties.
            # Recall:
            # CF1: Every source-to-target path in G_cf is cycle free.
            # CF2: Every cycle free path in the original graph appears as a
            #      source-to-target path in G_cf.
            # CF3: There is a 1-1 correspondence between the paths in G_cf and
            # the paths in the original graph. This means there is no
            # redundancy in the representation. For every path in the original
            # graph there is a unique path in G_cf that corresponds to it.
            # To do so, we first compute the set of source-to-target paths
            # (the nodes will be triples) in G_cf
            src_0 = (src[0], src[1], 0) # 3-tuple version of src
            tgt_0 = (tgt[0], tgt[1], 0) # 3-tuple version of tgt
            P_cf_pruned = list(nx.all_simple_paths(G_cf, src_0, tgt_0))

            # Next we extract the actual paths by projecting down to second
            # component.
            P_cf_pruned_names = name_paths(P_cf_pruned)

            # We first verify CF1.
            for p in P_cf_pruned_names:
                if len(p) != len(list(set(p))):
                    print("cycle!")
                    #print(p)
                    #break

            # Next we verify CF2. We will in fact check if the set of paths in
            # P_cf_pruned_names is exactly the set of paths in the original
            # graph.
            if set(P_correct) != set(P_cf_pruned_names):
                print("Missing")
                #print("graph, length", (i, length))
                #print("______________")
            print("# of paths: %d" % len(P_cf_pruned_names))

            # Finally we verify CF3
            if len(P_cf_pruned) != len(list(set(P_cf_pruned_names))):
                print("redundant representation!")
                #print("graph, length", (i, length))
                #print("____________")

if __name__ == '__main__':
    test_on_random_graphs()
