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
from .util import prune
import pickle
import numpy as np


def from_graph(g, source_name, target_name, path_length, fwd_reachset=None,
               back_reachset=None):
    pre_cfpg = pcf.from_graph(g, source_name, target_name, path_length,
                              fwd_reachset, back_reachset)
    return from_pre_cfpg(pre_cfpg, source_name, target_name, path_length)


def from_pre_cfpg(pre_cfpg, source_name, target_name, path_length):
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
    # Define old (2-tuple) and new (3-tuple) versions of src/tgt nodes
    src_2node = (0, source_name) # 2-tuple version of source
    src_3node = (0, source_name, 0) # 3-tuple version of source
    tgt_2node = (path_length, target_name) # 2-tuple version of target
    tgt_3node = (path_length, target_name, 0) # 3-tuple version of target
    # We first hardwire the contents of the dictionary for the level of the
    # target node: dic_CF[path_length]
    next_tgt = {tgt_3node: []}
    pred_tgt = {tgt_3node: pre_cfpg.graph.predecessors(tgt_2node)}
    t_cf_tgt = {tgt_3node: pre_cfpg.tags[tgt_2node]}
    dic_CF = {path_length: ([tgt_3node], next_tgt, pred_tgt, t_cf_tgt)}

    # If we were given an empty pre-CFPG, then the CFPG should also be empty
    if not pre_cfpg.graph:
        return nx.DiGraph()
    # Iterate from level n-1 (one "above" the target) back to the source
    for i in reversed(range(1, path_length)):
        # Get the information for level i+1 (one level closer to the target)
        V_ip1, next_ip1, pred_ip1, t_cf_ip1 = dic_CF[i+1]
        # Because we are working off of a non-empty pre-CFPG, we should never
        # end with a level in the graph with no nodes
        assert V_ip1 != []
        # TODO: Can V_current be replaced simply by the nodes in pre-CFPG at
        # level i?
        # TODO: Rename V_current -> V_i_old, V_i -> V_i_new?
        V_current = []
        for v in V_ip1:
            V_current.extend(pred_ip1[v])
            V_current = list(set(V_current))
        # V_current should never be empty by construction of the pre-CFPG
        assert V_current != []
        # Thus V_current is the set of nodes (which will be 2-tuples) at level
        # i to be processed. The converted  nodes (which will be 3-tuples,
        # including the copy number) will be binned into V_i.
        V_i, next_i, pred_i, t_cf_i = ([], {}, {}, {})
        # Now comes the heart of the construction. We take a node x in
        # V_current and split it into--in general--multiple copies to ensure
        # that if (u,v) is an edge in G_cf then the set of tags of u is
        # included in the set of tags of v
        for x in V_current:
            # X_ip1 is the set of nodes at the level i+1 to which x is connected
            # via the pred_ip1 function. These nodes, already processed, will
            # be 3-tuples. X_im1 are the set of predecessor nodes of x at the
            # level i-1. They are unprocessed 2-tuples.
            X_ip1 = [w for w in V_ip1 if x in pred_ip1[w]]
            X_im1 = pre_cfpg.graph.predecessors(x)
            assert X_ip1 != []
            # The actual splitting of node x and connect the resulting copies
            # of x to its neighbors above and below is carried out by the
            # _split_graph function, below.
            V_x, next_x, pred_x, t_cf_x = \
                    _split_graph(src_2node, tgt_2node, x, X_ip1, X_im1,
                                 t_cf_ip1, pre_cfpg)
            # We now extend V_i, next_i, pred_i and t_i in the obvious way.
            V_i.extend(V_x) # V_x contains the new, split versions of x
            next_i.update(next_x)
            pred_i.update(pred_x)
            t_cf_i.update(t_cf_x)
        dic_CF[i] = (V_i, next_i, pred_i, t_cf_i)
    # Finally we hardwire dic_CF[0]
    V_1 = dic_CF[1][0]
    V_0 = [src_3node]
    next_src = {src_3node: V_1}
    pred_src = {src_3node: []}
    t_cf_src = {src_3node: pre_cfpg.tags[src_2node]}
    dic_CF[0] = (V_0 , next_src, pred_src, t_cf_src)
    G_cf = _dic_to_graph(dic_CF)
    # Prune out possible unreachable nodes in G_cf
    nodes_prune = [v for v in G_cf
                     if (v != tgt_3node and G_cf.successors(v) == []) or
                        (v != src_3node and G_cf.predecessors(v) == [])]
    G_cf_pruned = prune(G_cf, nodes_prune, src_3node, tgt_3node)
    return CFPG(source_name, target_name, G_cf_pruned, path_length)


class CFPG(object):
    def __init__(self, source_name, target_name, graph, path_length):
        self.source_name = source_name
        self.source_node = (0, source_name, 0)
        self.target_name = target_name
        self.target_node = (path_length, target_name, 0)
        self.graph = graph
        self.path_length = path_length

    def enumerate_paths(self, names_only=True):
        paths = list(nx.all_simple_paths(self.graph, self.source_node,
                                         self.target_node))
        if names_only:
            return self._name_paths(paths)
        else:
            return paths

    @staticmethod
    def _name_paths(paths):
        return [tuple([node[1] for node in path]) for path in paths]


def _split_graph(src, tgt, x,  X_ip1, X_im1, t_cf, pre_cfpg):
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
        X_wx =  set(t_cf[w]) & set(pre_cfpg.tags[x])
        N_wx = list(X_wx)
        # TODO: Reimplement pruning so as to avoid inducing a subgraph?
        g_wx = pre_cfpg.graph.subgraph(N_wx)
        nodes_prune = [v for v in g_wx
                         if (v != x and g_wx.successors(v) == []) or
                            (v != src and g_wx.predecessors(v) == [])]
        g_wx_pruned = prune(g_wx, nodes_prune, src, x)
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


