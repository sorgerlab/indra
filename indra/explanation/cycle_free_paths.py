import random
import itertools
from copy import copy
import networkx as nx
from .paths_graph import get_edges, get_reachable_sets, paths_graph

"""
As we know, a paths-graph, say G_n, has the property that every node in G_n
lies on a path of length n from source to target. These paths may not be
cycle-free. Our goal is to transform G_n,  into a new graph G_cf such that in
G_cf every node lies on a cycle-free path from the source to a target. However
it will NOT be the case -despite the misleading name G_cf- that G_cf contains
only cycle-free paths. However we will be able to "easily" sample cycle-free
paths from G_cf without any backtracking. My hope is that without too large a
blow-up, G_cf can be transformed to G_CF which contains only cycle-free paths.
In fact I have a rough idea for doing this but haven't tested it out.

In the present experiment we work with the korkut_im.sif example. We set n = 10
(also equal to max_depth in the get_reachable_sets function).  We will go from
G_10( called pg_raw below) to G_cf in stages. At stage i we process all the
nodes at level i.  We start from G_0  at level 0. The first step is  special in
which obtain G_1 by processing the nodes (0, source) and  (10, target). Then by
processing the the nodes of G_1 at level 1 we will obtain G_2 etc. At the end
of this process we will set G_cf = G_8

The input to the following forward reachset computation will be the graph G_i
at the i-th stage for 1 <= j <= 8 together with one of its nodes at level j.
This forward reach set will be pruned at the next step as explained below.
"""

def forward(v, src, tgt, H):
    j = v[0]
    L = {}
    L[j] = [v]

    h = nx.DiGraph()

    for k in range(j+1, 10):
        for v in L[k-1]:
            h.add_edges_from(H.out_edges(v))
        L[k] = [w for w in h if w[0] == k]
    #h.add_nodes_from(V)
    return h

"""
We also compute the backward version.  We will then attach the pruned forward
reach set to the backward reachset to obtain G_(j+1, v). By patching together
all the graphs in {G_{j+1, v}}_{v in j-th level of G_j} we will obtain G_j.
Remark: The actual names in the code are slightly different. This should be
fixed. Sorry!
"""

def backward(v, src, tgt, H):
    j = v[0]
    L = {}
    L[j] = [v]
    J =  list(reversed(range(0, j)))
    h = nx.DiGraph()
    for k in J:
        for v in L[k+1]:
            h.add_edges_from(H.in_edges(v))
        L[k] = [w for w in h if w[0] == k]
    return h


"""
forward(v, source, target, G_j) will return the subgraph of G_j, say G_j_v_f,
that is forward reachable from v in G_j.

We prune G_j_v_f through the following iterative procedure.  Let X_j_v be the
set of nodes in G_j_v_f -other than v- that have the same name as v.  In other
words u[1] = v[1] for every u in X_j_v.

The function span(U, src, tgt, g) takes as inputs (a) g the current (possibly
partially pruned) graph (b) U, the set of nodes U that can be pruned.

(i) Remove all the nodes in U.

(ii) Remove all the in_coming and out_going edges (relative to G_j_v_f)  of
these nodes.

(iii) In the resulting graph, identify Y, the set of nodes  have lost all their
in_coming edges or out_going edges.

Return g_prune, the partially pruned graph (obatined through steps (i) and
(ii)) and Y, the candidiate set of nodes for the next iteration of pruning.

prune(Y, src, tgt, g) applies span() repeatedly till there are no more nodes to
be pruned. We call prune() with Y = X_j_v and g = G_j_v_f

In the graph resulting from this pruning process, suppose both v and and tgt
(ie. the node (9, target)) have survived. Since we will be going from lower
levels to higher levels (yes, I can now think forwards!) we are guaranteed to
be able to trace a cycle-free path from (0, source) to v. Further, thanks to
the pruning we are guaranteed to be able to go from v to tgt without hitting a
cycle involving v[1] in the current stage.

If v or tgt gets eliminated then we can conclude that there are no cycle free
paths passing through v in G_j and proceed accordingly.
"""

"""
def span(U, g):
    g_current = nx.DiGraph()
    g_current.add_edges_from(g.in_edges(U))
    g_current.add_edges_from(g.out_edges(U))

    g.remove_edges_from(g_current.edges())
    h = nx.DiGraph()
    h.add_edges_from(g.edges())

    Y = []
    L = [n for n in g_current.nodes_iter()]
    E = [e for e in g_current.edges_iter()]
    for v in L:
        if v in h.nodes():
            I_v = set(h.in_edges(v)) 
            O_v = set(h.out_edges(v))
            if (I_v <= set(E)) or (O_v <= set(E)):
                Y.append(v)
    return (h, Y)


def prune(U, g):
    while U != []:
        g_prune, Y = span(U, g)
        g = g_prune.copy()
        U = Y
    return g_prune
"""

def prune(pg, nodes_to_prune, path_length):
    # Make a copy of the graph
    pg_pruned = pg.copy()
    while nodes_to_prune:
        pg_pruned.remove_nodes_from(nodes_to_prune)
        # Make a list of nodes whose in or out degree is now 0 (making
        # sure to exclude the source and target, whose depths are at 0 and
        # path_length, respectively)
        no_in_edges = [node for node, in_deg in pg_pruned.in_degree_iter()
                        if in_deg == 0 and node[0] != 0]
        no_out_edges = [node for node, out_deg in pg_pruned.out_degree_iter()
                        if out_deg == 0 and node[0] != path_length]
        nodes_to_prune = set(no_in_edges + no_out_edges)
    return pg_pruned

"""
For the purposes of sampling cycle-free paths we attach a set of tags to each
node when we compute G_j+1 from G_j.  Basically we add the tag v[1] to each
node in G_j+1 that can be reached from v in the forward direction. This tag
tells that in tracing path starting from v as long as we stick to nodes whose
tag set contains v, we will be ok. 

pg_raw is our path graph. As mentioned above the first step is special. We
eliminate all cycles involving src as well as all cycles involving tgt. To
prime the sampling procedure we add the tag 'source' to every node in the
pruned graph; except to src whose tag set will be [] """

def PG_0(src, tgt, pg_raw):
    g = pg_raw.copy()
    import ipdb; ipdb.set_trace()
    """ while debugging I got into trouble by not keeping track of different
    versions of the same graph. Hence playing it safe above """ 
    """ First identify the nodes to be pruned. They are just nodes whose names
    are either 'source' or 'target' """
    Z_raw = []
    for v in g.nodes_iter():
        if (v != src) & (v != tgt) & ((v[1] == src[1]) or (v[1] == tgt[1])):
            Z_raw.append(v)
    Z_raw = list(set(Z_raw))

    """ If there are no nodes to be pruned (this is the case in our current
    example) just add the tag [source] to every node other than src. The node
    src gets the tag []."""
    if Z_raw == []:
        # Dictionary for storing tags, indexed by node. Tags consist of a
        # list of node names (without depths)
        tags_0 = {}
        for v in pg_raw.nodes_iter():
            if v == src:
                tags_0[v] = []
            else:
                tags_0[v] = [src[1]]
        # Because we didn't prune out any nodes, we can return the original
        # paths graph
        return (pg_raw, tags_0)
    else:
        g_pruned  = prune(Z_raw, src, tgt, g)
        # If the source or target gets pruned then there are no cycle free
        # paths. Hence we return the degenerate (graph, tags) pair and
        # propagate it through the remaining stages.
        if (src in g_pruned) or (tgt in g_pruned):
            g_empty = nx.DiGraph()
            tags_empty = {}
            return (g_empty, tags_empty)
        else:
            tags_0 = {}
            for v in g_pruned.nodes_iter():
                if v == src:
                    tags_0[v] = []
                else:
                    tags_0[v] = [src[1]]
        return (g_pruned, tags_0)

""" Finally we compute each G_j for 1 <= j <= 8 together with tag sets """

def PG(src, tgt, pg_0):
    dic_PG = {}
    dic_PG[0] = pg_0
    for k in range(1, 9):
        H = dic_PG[k-1][0].copy()
        tags = dic_PG[k-1][1].copy()
        """ Check if we have already detected there are no cycle free paths. If
        so just propagate this information """
        if H.nodes() == []:
            dic_PG[k] = dic_PG[k-1]
        else:
            """ Identify the nodes at level k in G_(k-1) """
            X = [v for v in H.nodes_iter() if v[0] == k]
            """ we will track the (g_x, tags_x) pairs constributed by each x
            through dic_X """
            dic_X = {}
            for x in X:
                tags_x = {}
                g_x_b = backward(x, src, tgt, H)
                g_x_f = forward(x, src, tgt, H)
                g_x = nx.DiGraph()
                g_x.add_edges_from(g_x_b.edges())
                g_x.add_edges_from(g_x_f.edges())
                Z_x = [v for v in g_x_f if ((v[1] == x[1]) & (v[0] != k))]
                """  If Z_x is null then just add the tag 'x' all the nodes in
                g_x_f but not to x"""
                if Z_x == []:
                    for v in g_x.nodes_iter():
                        if v[0] > k:
                            D = tags[v]
                            D.append(x[1])
                            tags_x[v] = D
                        else:
                            tags_x[v] = tags[v]
                    dic_X[x] = (g_x, tags_x)
                else:
                    """ carry out the pruning """
                    g_x_prune = prune(Z_x, src, tgt, g_x)
                    """ if tgt or x gets pruned then x will contribute nothing
                    to G_k """
                    if (tgt not in g_x_prune) or (x not in g_x_prune):
                        pass
                    else:
                        """ otherwise add the tag x to the nodes in the strict
                        future of x. update dic_X """
                        for v in g_x_prune.nodes_iter():
                            if v[0] > k:
                                D = tags[v]
                                D.append(x[1])
                                tags_x[v] = D
                            else:
                                tags_x[v] = tags[v]
                    dic_X[x] = (g_x_prune, tags_x)
            """ we can now piece together the pairs in dic_X  to obatin (G_k,
            tags_k) """
            H_k = nx.DiGraph()
            tags_k = {}
            for x in X:
                h_x = dic_X[x][0]
                H_k.add_edges_from(h_x.edges())
            for v in H_k.nodes_iter():
                t = []
                for x in X:
                    if v in dic_X[x][0]:
                        tags_x = dic_X[x][1]
                        t.extend(tags_x[v])
                t = list(set(t))
                tags_k[v] = t
        dic_PG[k] = (H_k, tags_k)
    return dic_PG

""" The sampling procedure simply uses the tag sets to trace out cycle-free
paths. The basic idea is if we have reached a node v via the path p then we can
choose the successor u of v as the next node only if p appears in the tag set
of u. """

def cf_succ(H,t, path, v):
    succ = []
    for u in H.successors(v):
        if set(path) <= set(t[u]):
            succ.append(u)
    w = random.choice(succ)
    return w

def cf_sample_single_path(src, tgt, H,t):   
        path = [src[1]]
        current = src
        while current != tgt:
            next = cf_succ(H,t, path, current)
            """ a sanity check; since I have not stree-tested the code yet """
            if next[1] in path:
                print(path.append(next))
                print("error in the computation of G_cf")
                break
            else:
                path.append(next[1])
                current = next
        return tuple(path)

def cf_sample_many_paths(src,tgt, H, t, n):
    P = []
    for i in range(0, n):
        p = cf_sample_single_path(src, tgt, H,t)
        P.append(p)
    return P


if __name__ == '__main__':
    G_0 = get_edges('korkut_im.sif')
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'

    (f_level, b_level)  =  get_reachable_sets(G_0, source, target,
                                              max_depth=10, signed=False)
    length = 9

    pg_raw = paths_graph(G_0, source, target, length, f_level, b_level,
                    signed=False, target_polarity=0)

    src = (0, source)
    tgt = (9, target)
    pg_0 = PG_0(src,tgt,pg_raw)

    dic_PG = PG(src,tgt,pg_0)
    G_cf, T = dic_PG[8]
    P = cf_sample_many_paths(src,tgt,G_cf, T, 1000)
    print(len(list(set(P))))
    #print("--- %s seconds ---" % (time.time() - start_time))
