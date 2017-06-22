# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:28:25 2017

@author: thiagu

junk
"""

import random
import itertools
import networkx as nx

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

""" We can now construct for any chosen length l , a graph that contains all paths of length l from 
source to target. l is the number of edges encountered onhe path. We assume 2 <= l <= depth. The idea is to sample for paths using this graph. Once a path has been sampled we can prune cycles from it.
We can also inject polrities onto the edges of the path and thus compute the parities of the nodes encountered on the path. """ 

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

def pre_7(v):
    X = []
    for u in N:
        if (u,v) in A:
            X.append(u)
    return X

def sample_next(current):
    X = pre_7(current)
    c = random.randint(0,len(X) - 1)
    new = X[c]
    return new

def sample_path(g, source, target, length, f_level, b_level):
    path = []
    path.append((0, source))
    while len(path) < length:
        current_new = path[-1]
        u = sample_next(current_new)
        path.append(u)
    path.append((length, target))
    path_pruned = []
    for i in range(len(path)):
        path_pruned.append(path[i][1])
    return path_pruned


if __name__ == '__main__':
    g = get_edges('korkut_im.sif')
    #g = get_edges('test_graph.sif')
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'
    #source = 'A'
    #target = 'D'
    target_polarity = 1

    # fix the maximuma dpth to which we want to explore backwards from the
    # traget and forwaards from the source
    print("Computing forward and backward reach sets...")
    max_depth = 10

    # Compute backward level sets
    b_level = {}
    b_level[0] = set([(target, 0)])

    for i in range(1, max_depth):
        X = []
        for v, node_polarity in b_level[i-1]:
            preds = g.predecessors(v)
            for u in preds:
                edge_polarity = g.get_edge_data(u, v)['polarity']
                cum_polarity = (node_polarity + edge_polarity) % 2
                X.append((u, cum_polarity))
        if not X:
            break
        b_level[i] = set(X)

    # Compute forward level sets
    f_level = {}
    f_level[0] = set([(source, 0)])

    for i in range(1, max_depth):
        X = []
        for u, node_polarity in f_level[i-1]:
            succs = g.successors(u)
            for v in succs:
                edge_polarity = g.get_edge_data(u, v)['polarity']
                cum_polarity = (node_polarity + edge_polarity) % 2
                X.append((v, cum_polarity))
        if not X:
            break
        f_level[i] = set(X)

    # Compute path graph for a specific path length
    length = 10
    print("Computing path graph of length %d" % length)
    pg = paths_graph(g, source, target, target_polarity, length, f_level,
                     b_level)
    print("Drawing Graphviz graph")
    ag = nx.nx_agraph.to_agraph(pg)
    ag.draw('gop.pdf', prog='dot')
    # -----
    """
    Paths = []
    while len(Paths) < 10:
        p = sample_path(g, source, target, length, f_level, b_level)
        Paths.append(p)


    g = open('Paths_7.txt', 'w')

    for i in range(10):
        print >> g, 'Paths_7[',i,'] =' 
        for j in range(len(Paths[i])): 
            print >> g, Paths[i][j]
        print >> g, ""
    g.close()
    """
