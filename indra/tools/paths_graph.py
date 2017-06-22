import random
import itertools
import logging
import networkx as nx

logger = logging.getLogger('paths_graph')

def get_edges(sif_file):
    edges = []
    with open(sif_file, 'rt') as f:
        for line in f.readlines():
            u, polarity, v = line.strip().split(' ')
            if u == v:
                pass
            else:
                edges.append((u, v, {'polarity': int(polarity)}))
    g = nx.DiGraph()
    g.add_edges_from(edges)
    return g

""" We can now construct for any chosen length l , a graph that contains all
paths of length l from source to target. l is the number of edges encountered
onhe path. We assume 2 <= l <= depth. The idea is to sample for paths using
this graph. Once a path has been sampled we can prune cycles from it.  We can
also inject polrities onto the edges of the path and thus compute the parities
of the nodes encountered on the path. """ 

def paths_graph(g, source, target, target_polarity, length, f_level, b_level):
    level = {}
    level[0] = set([(target, target_polarity)])
    level[length] = set([(source, 0)])
    # If the target polarity (from source to target) is odd, flip the polarity
    # of all nodes in the graph
    if target_polarity == 1:
        b_level_polar = {}
        for i in range(0, len(b_level)):
            polar_set = set()
            for (u, w) in b_level[i]:
                w_flipped = (w + 1) % 2
                polar_set.add((u, w_flipped))
            b_level_polar[i] = polar_set
    else:
        b_level_polar = b_level

    for i in range(1, length):
        b_reach_set = b_level_polar[i]
        f_reach_set = f_level[length - i]
        Z = set(f_reach_set) & set(b_reach_set)
        level[i] = Z
    V_gop = {}
    for i in range(0,length+1):
        V_gop[i] = list(itertools.product([i], level[i])) 
    E_gop = set()
    for i in range(0, length):
        X = set() # list of edges at this level
        for u, v in itertools.product(V_gop[i+1], V_gop[i]):
            u_name, u_pol = u[1]
            v_name, v_pol = v[1]
            if (u_name, v_name) in g.edges():
                edge_polarity = g.get_edge_data(u_name, v_name)['polarity']
                # Look for an edge the flips or doesn't flip the polarity
                # of the path depending on what we see in the cumulative
                # polarities
                if (u_pol == v_pol and edge_polarity == 0) or \
                   (u_pol != v_pol and edge_polarity == 1):
                    X.add((u, v))
        E_gop |= X
    path_graph = nx.DiGraph()
    path_graph.add_edges_from(E_gop)
    return path_graph


def sample_path(g, source, target, length):
    path = []
    source_node = (length, (source, 0))
    path.append(source_node)
    assert source_node in g
    for i in range(length):
        current_node = path[-1]
        succs = g.successors(current_node)
        v = random.choice(succs)
        path.append(v)
    return path


def get_reach_sets(g, max_depth=10):
    logger.info("Computing forward and back reach sets")
    # Backward level sets
    b_level = {}
    b_level[0] = set([(target, 0)])
    # Forward level sets
    f_level = {}
    f_level[0] = set([(source, 0)])

    directions = ((f_level, lambda u: [((u, v), v) for v in g.successors(u)]),
                  (b_level, lambda v: [((u, v), u) for u in g.predecessors(v)]))
    for level, edge_func in directions:
        for i in range(1, max_depth):
            reachable_set = set()
            for node, node_polarity in level[i-1]:
                for (u, v), reachable_node in edge_func(node):
                    edge_polarity = g.get_edge_data(u, v)['polarity']
                    cum_polarity = (node_polarity + edge_polarity) % 2
                    reachable_set.add((reachable_node, cum_polarity))
            if not reachable_set:
                break
            level[i] = reachable_set
    return (f_level, b_level)


if __name__ == '__main__':
    g = get_edges('korkut_im.sif')
    #g = get_edges('test_graph.sif')
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'
    #source = 'A'
    #target = 'D'
    target_polarity = 0

    # fix the maximuma dpth to which we want to explore backwards from the
    # traget and forwaards from the source
    print("Computing forward and backward reach sets...")
    max_depth = 10
    f_level, b_level = get_reach_sets(g, max_depth)

    # Compute path graph for a specific path length
    length = 6
    print("Computing path graph of length %d" % length)
    pg = paths_graph(g, source, target, target_polarity, length, f_level,
                     b_level)
    print("Drawing Graphviz graph")
    ag = nx.nx_agraph.to_agraph(pg)
    ag.draw('gop.pdf', prog='dot')

    paths = []
    while len(paths) < 10:
        p = sample_path(pg, source, target, length)
        paths.append(p)

    print(paths)

