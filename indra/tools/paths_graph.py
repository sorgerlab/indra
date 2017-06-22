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
            u, weight, v = line.strip().split(' ')
            if u == v:
                pass
            else:
                edges.append((u, v, {'weight': int(weight)}))
    g = nx.DiGraph()
    g.add_edges_from(edges)
    return g

""" We can now construct for any chosen length l , a graph that contains all paths of length l from 
source to target. l is the number of edges encountered onhe path. We assume 2 <= l <= depth. The idea is to sample for paths using this graph. Once a path has been sampled we can prune cycles from it.
We can also inject polrities onto the edges of the path and thus compute the parities of the nodes encountered on the path. """ 

def GoP(length, f_level, b_level, source, target):
    level = {}
    level[0] = [target]
    level[length] = [source]
    for i in range(1,length):
        X = b_level[i]
        Y = f_level[length - i]
        Z = set(X) & set(Y)
        level[i] = Z
    V_gop = {}
    for i in range(0,length+1):
        V_gop[i] = list(itertools.product([i], level[i])) 
    E_gop = {}
    for i in range(0, length):
        X = []
        for e in itertools.product(V_gop[i+1], V_gop[i]):
            if (e[0][1], e[1][1]) in E:
                X.append(e)
        E_gop[i] = X

    V_out = []
    for i in range(0, length + 1):
        V_out.extend(V_gop[i])
    E_out = []
    for i in range(0, length):
        E_out.extend(E_gop[i])   

    return level, V_out, E_out


def GoP_alt(g, source, target, length, f_level, b_level):
    level = {}
    level[0] = [target]
    level[length] = [source]
    for i in range(1,length):
        X = b_level[i]
        Y = f_level[length - i]
        Z = set(X) & set(Y)
        level[i] = Z
    V_gop = {}
    for i in range(0,length+1):
        V_gop[i] = list(itertools.product([i], level[i])) 
    E_gop = set()
    for i in range(0, length):
        X = set() # list of edges at this level
        for u, v in itertools.product(V_gop[i+1], V_gop[i]):
            if (u[1], v[1]) in g.edges():
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
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'

    # fix the maximuma dpth to which we want to explore backwards from the
    # traget and forwaards from the source
    max_depth = 10

    # Compute backward level sets
    b_level = {}
    b_level[0] = [target]

    for i in range(1, max_depth):
        X = []
        for v in b_level[i-1]:
            X.extend(g.predecessors(v))
        b_level[i] = set(X)

    # Compute forward level sets
    f_level = {}
    f_level[0] = [source]

    for i in range(1, max_depth):
        X = []
        for v in f_level[i-1]:
            X.extend(g.successors(v))
        f_level[i] = list(set(X))

    # Compute path graph for a specific path length
    length = 10
    gop_n = GoP_alt(g, source, target, length, f_level, b_level)
    #ag = nx.nx_agraph.to_agraph(gop_n)
    #ag.draw('gop.pdf', prog='dot')
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
