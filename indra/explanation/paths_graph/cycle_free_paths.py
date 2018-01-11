"""
We construct a representation of cycle_free paths of a fixed length. This fixed
length will often not be mentioned in what follows.  We call our representation
"the cycle_free paths graph". Below it is the graph G_cf (actually it is
G_cf_pruned but for now it will be convenient to ignore this distinction).

G_cf is required to have three properties.

* CF1: Every source-to-target path in G_cf is cycle free.
* CF2: Every cycle free path in the original graph appears as a
  source-to-target path in G_cf.
* CF3: There is a 1-1 correspondence between the paths in G_cf and the paths in
  the original graph. This means there is no redundancy in the representation.
  For every path in the original graph there is a unique path in G_cf that
  corresponds to it.

These 3 conditions will ensure that we can sample paths in the original graph
faithfully by sampling paths in G_cf. We can also perform graph theoretic
operations on G_cf to simulate useful operations on the set of paths in the
original graph.

The starting point is the paths graph (pg_raw below) that represents "all"
paths (cycle free or not) of the given fixed length from source to target.

Then using an initial iterative procedure we prune away junk nodes (that cannot
appear on any cycle free path from source to target) and more importantly tag
each node with its cycle free history. More precisely if u is in tags[v] then
we are guaranteed that every path from u to v that involves only nodes
appearing in tags[v] will be v-cycle_free. In other words the name of v (i.e.
v[1]) will not appear in the path. Further, it will also be u-cycle free.  Note
however tags[u] may contain a node that has the same name as that of v.  Indeed
this the crux of the problem.

Moving on, this tagged path graph is named G_0 and the associated tags map is
named T_0 below.

G_cf is computed by refining G_0. But first let us consider why G_0 is not an
ideal representation of the set of cycle free paths of a fixed length.  First,
G_0 does not have the property (CF1) (though it does have the properties (CF2)
and (CF3)). As a result one can't just walk through the graph from source to
node and generate a cycle free path. Instead one must use a sampling method
with memory to generate cycle free paths. In particular if one has reached the
node u via the path p and v is a successors of u then one can extend p by
moving to v only if p is contained in T_0[v]. Thus whether the move along the
edge (u,v) is conditioned by the memory of how u was reached. Further, one can
get stuck while using this sampling procedure. Hence it is not clear whether
one is sampling the set of paths of interest in a faithful fashion. More
importantly it is not clear how one can perform graph theoretic operations on
G_0 to simulate operations on the set of cycle fre paths of interest. We will
however keep in mind that G_0 together with its path sampling procedure is  a
useful tool to have around.

Constructing G_cf by refining G_0 may be viewed as synthesizing a memoryless
strategy for generating cycle free paths. In other words, if (u,v) is an edge
in G_cf then no matter how we have reached u we must be able to transition to
v. A necessary condition that will enable this is to ensure that the set of
tags of u (T_cf[u]) is included in the set of tags of v (T_cf[v]) in G_cf. The
challenge is to achieve this while ensuring that the properies (CF1), (CF2) and
(CF3) are met.

We explain below the detailed construction of G_cf from this perpective.
"""

import random
import itertools
import networkx as nx
from . import paths_graph
from . import pre_cfpg as pcf
import pickle
import numpy as np


def PG_cf(src, tgt, g, t):
    """Generate a cycle free paths graph (CFPG).

    Implements the major step (the outer loop) for constructing G_cf. We do so
    by computing dic_CF, a dictionary based version of G_cf.  dic_CF[i] will be
    a quadruple of the form (V_i, next_i, pred_i, t_i).

    V_i will be the set of nodes at level i.

    A node--after dic_CF[i] has been computed--will be of the form (i, n, c)
    where i is the level, n is the name and c is the copy number of the node
    (i,n) in G_0. In other words, each node in G_0 will be split into one or
    more copies to implement our memoryless sampling strategy.

    next_i is the successor relation for the CFPG.

    pred_i[v] is the set of predecessors of v in V_i. The construction
    proceeds from the target to source. At stage i of the construction we
    convert nodes of the form (i, n) into nodes of the form (i,n,c). For any
    such new node pred_i[v] will be nodes of the form (i-1,n) at level i-1.

    t_i[v] will be the new tags of the node v. They will be pairs of the form
    (j,n). In other words their type will be the same as of T_0.  (Note: In
    T_0, I assign nodes of G_0 as tags rather than their names. This turns out
    to be convenient for the construction of G_cf)

    Once the construction of PG_cf is complete we will no onger 
    require pred_i and t_i.
    """
    # We first hardwire dic_CF[tgt[0]].
    tgt_0 = (tgt[0], tgt[1], 0) # 3-tuple version of tgt
    next_tgt = {tgt_0: []}
    pred_tgt = {tgt_0: g.predecessors(tgt)}
    t_cf_tgt = {tgt_0: t[tgt]}
    dic_CF = {tgt[0]: ([tgt_0], next_tgt, pred_tgt, t_cf_tgt)}

    # If we were given an empty pre-CFPG, then the CFPG should also be empty
    if not g:
        return nx.DiGraph()
    # Iterate from level n-1 (one "above" the target) back to the source
    for i in reversed(range(1, tgt[0])):
        # Get the information for level i+1 (one level closer to the target)
        V_ip1, next_ip1, pred_ip1, t_cf_ip1 = dic_CF[i+1]
        # Because we are working off of a non-empty G_0, we should never end
        # with a level in the graph with no nodes
        assert V_ip1 != []
        # TODO: Can V_current be replaced simply by the nodes in g at level i?
        # TODO: Rename V_current -> V_i_old, V_i -> V_i_new?
        V_current = []
        for v in V_ip1:
            V_current.extend(pred_ip1[v])
            V_current = list(set(V_current))
        # V_current should never be empty by construction of G_0
        assert V_current != []
        # Thus V_current is the set of nodes (which will be 2-tuples) at level
        # i to be processed. The converted  nodes(which will be 3-tuples,
        # including the copy number) will be binned into V_i.
        V_i = []
        next_i = {}
        pred_i = {}
        t_cf_i = {}
        # Now comes the heart of the construction. We take a node x in
        # V_current and split it into -in general- multiple copies to ensure
        # that if (u,v) is an edge in G_cf then the set of tags of u is
        # included in the set of tags of v
        for x in V_current:
            # X_ip1 is the set of nodes at the level i+1 to which x is connected
            # via the pred_ip1 function. These nodes, already processed, will
            # be 3-tuples. X_im1 is the nodes at the level i-1. They are
            # unprocessed 2-tuples.
            X_ip1 = [w for w in V_ip1 if x in pred_ip1[w]]
            X_im1 = g.predecessors(x)
            assert X_ip1 != []
            # The actual splitting up and how to connect the resulting copies
            # of x to its neighbors above and below is carried out by the
            # _split_graph function, below.
            V_x, next_x, pred_x, t_cf_x = \
                    _split_graph(src, tgt, x,  X_ip1, X_im1, t_cf_ip1, t, g)
            # We now extend V_i, next_i, pred_i and t_i in the obvious way.
            V_i.extend(V_x)
            next_i.update(next_x)
            pred_i.update(pred_x)
            t_cf_i.update(t_cf_x)
        dic_CF[i] = (V_i, next_i, pred_i, t_cf_i)
    # Finally we hardwire dic_CF[0]
    V_1 = dic_CF[1][0]
    src_0 = (src[0], src[1], 0) # 3-tuple version of src
    V_0 = [src_0]
    next_src = {src_0: V_1}
    pred_src = {src_0: []}
    t_cf_src = {src_0: t[src]}
    dic_CF[0] = (V_0 , next_src, pred_src, t_cf_src)
    G_cf = _dic_to_graph(dic_CF)
    # Prune out possible unreachable nodes in G_cf
    nodes_prune = [v for v in G_cf
                     if (v != tgt_0 and G_cf.successors(v) == []) or
                        (v != src_0 and G_cf.predecessors(v) == [])]
    G_cf_pruned = pcf._prune(G_cf, nodes_prune, src_0, tgt_0)
    return G_cf_pruned


def _split_graph(src, tgt, x,  X_ip1, X_im1, t_cf, t, g):
    """Splits a node x from G_0 into multiple copies for the CFPG.

    The nodes in X_ip1 represent the possible successor nodes to x in the CFPG.
    For each successor w of x in X_ip1, we first obtain the set of possible
    antecedent nodes lying on paths from the source up to the edge x->w. We
    obtain this by finding the intersection between the tags of x and the tags
    of w. This is the set X_wx below for each w in X_ip1.

    However X_wx is the set of nodes (in G_0) from which we can reach x->w
    without encountering x[1] AND without encountering w[1]. As a result some
    nodes in X_wx may be isolated.  Hence we prune them away.
    """
    V_x = []
    next_x = {}
    pred_x = {}
    t_x = {}
    S_ip1 = {}
    for w in X_ip1:
        X_wx =  set(t_cf[w]) & set(t[x])
        N_wx = list(X_wx)
        # TODO: Reimplement pruning so as to avoid inducing a subgraph?
        g_wx = g.subgraph(N_wx)
        nodes_prune = [v for v in g_wx
                         if (v != x and g_wx.successors(v) == []) or
                            (v!= src and g_wx.predecessors(v) == [])]
        g_wx_pruned = pcf._prune(g_wx, nodes_prune, src, x)
        # If the pruned graph still contains both src and x itself, there is
        # at least one path from the source to x->w. The nodes in this subgraph
        # constitute the new set of tags of the copy of x that lies on a path
        # between src and w.
        if x in g_wx_pruned and src in g_wx_pruned:
            s = frozenset(g_wx_pruned.nodes())
            S_ip1[w] = s
    S = set(S_ip1.values())
    # Each element of the set S will be a unique, (frozen) set of tags. We will
    # create one copy x_r of x for each unique tag set r in S, and we assign r
    # to be the set of tags of the new, split node x_r. The successors of x_r
    # are assembled using S_ip1; pred is defined in the expected way using
    # X_im1.
    for c, r in enumerate(S):
        x_c = (x[0], x[1], c)
        V_x.append(x_c)
        next_x[x_c] = [w for w in S_ip1.keys() if r == S_ip1[w]]
        pred_x[x_c] = [u for u in X_im1 if u in r]
        t_x[x_c] = r
    return (V_x, next_x, pred_x, t_x)


def _dic_to_graph(dic):
    G = nx.DiGraph()
    E = []
    for k in dic.keys():
        V_k = dic[k][0]
        next_k = dic[k][1]
        for v in V_k:
            E_v = list(itertools.product([v], next_k[v]))
            E.extend(E_v)
    G.add_edges_from(E)
    return G


def name_paths(Q):
    Q_names = []
    for q in Q:
        q_names = []
        for j in range(len(q)):
            q_names.append(q[j][1])
        Q_names.append(tuple(q_names))
    return Q_names


def sample_single_path(src_0, tgt_0, dic):
    p = [src_0]
    current = src_0
    while current != tgt_0:
        nxt = dic[current[0]][1]
        succ = random.choice(nxt[current])

        p.append(succ)
        current = succ
    return tuple(p)


def sample_many_paths(src_0, tgt_0, dic, n):
    P = []
    for i in range(0, n):
        p = sample_single_path(src_0, tgt_0, dic)
        P.append(p)
    return P


if __name__ == '__main__':
    # We use 25 randomly generated graphs for testing the algorithm
    dic_rg_long = pickle.load(open("save.dic_rg_long", "rb"))
    #for i in range(25):
    min_depth = 5
    max_depth = 5
    for i in range(1):
        G_i, source, target = dic_rg_long[i]
        print("graph#, no of nodes, no of edges",
              i, len(G_i), len(G_i.edges()))
        (f_level, b_level)  = \
                paths_graph.get_reachable_sets(G_i, source, target,
                        max_depth=max_depth, signed=False)
        # Try different path lengths
        for length in range(min_depth, max_depth+1):
            print("paths of length", length)
            print("_______")
            # For validation, we compute explicitly the set of paths in the
            # original graph of a fixed length
            P = list(nx.all_simple_paths(G_i, source, target, length+1))
            # Filter to paths of this length
            P_correct = [tuple(p) for p in P if len(p) == length+1]
            # Generate the raw paths graph
            pg_raw = paths_graph.paths_graph(G_i, source, target, length,
                                             f_level, b_level, signed=False,
                                             target_polarity=0)
            src = (0, source)
            tgt = (length, target)
            # The "pre" CFPG
            dic_PG = pcf.from_pg(pg_raw, src, tgt, length)
            # FIXME: Why isn't the target in dic_PG[len-1]?
            G_0, T_0 = dic_PG[len(dic_PG)-1]
            # The above is the output of the iterative method
            T_0[tgt].append(tgt)
            G_cf = PG_cf(src,tgt,G_0,T_0)

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
                    print(p)
                    break

            # Next we verify CF2. We will in fact check if the set of paths in
            # P_cf_pruned_names is exactly the set of paths in the original
            # graph.
            if set(P_correct) != set(P_cf_pruned_names):
                print("Missing")
                print("graph, length", (i, length))
                print("______________")
            print("# of paths: %d" % len(P_cf_pruned_names))

            # Finally we verify CF3
            if len(P_cf_pruned) != len(list(set(P_cf_pruned_names))):
                print("redundant representation!")
                print("graph, length", (i, length))
                print("____________")
