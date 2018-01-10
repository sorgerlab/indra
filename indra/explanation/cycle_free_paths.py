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
from copy import copy, deepcopy
import networkx as nx
from explanation import paths_graph
import cPickle as pickle


def cycle_free_paths_graph(pg, source, target, path_length):
    """Compute a pre cycle free paths graph.

    Starting from the "raw" (i.e., containing cycles) paths graph, and given a
    target path length n, the algorithm iterates over each "level" in the graph
    0 <= k <= n where level 0 consists only of the source node and level n
    consists only of the target.

    Each level k consists of a set of nodes, X; we examine each node x in X and
    identify the subset of nodes that are reachable in both the forward and
    backward directions from x. If any of the nodes in the forward reach
    subgraph contain x itself (but at a different depth), this represents a
    cyclic path back through x that is then pruned.

    Each node x therefore defines its own subgraph of cycle free paths, g_x.
    After iterating over all x in X, we combine these subgraphs into the
    (in-progress) cycle free paths graph H_k. H_k therefore consists of the
    superset of nodes of all the subgraphs g_x for level k. When merging these
    subgraphs we prevent the re-introduction of cyclic paths by annotating each
    node in the graph with a list of "tags". The tags for any given node
    consist of a list of nodes lying at prior (upstream) levels. Therefore
    during sampling, transitions from an upstream node to a downstream node are
    only permissible if all nodes in the path up to a certain level are
    contained in the tag set of the downstream node.

    Parameters
    ----------
    pg : networkx.DiGraph()
        "Raw" (contains cycles) paths graph as created by
        :py:func:`indra.explanation.paths_graph.paths_graph`.
    source : tuple
        Source node, of the form (0, source_name).
    target : tuple
        Target node, of the form (target_depth, source_name).
    path_length : int
        Desired path length.

    Returns
    -------
    tuple : (networkx.DiGraph(), dict)
        The initialized, not-yet cycle free paths graph consists of the
        paths graph remaining after cycles through the source or target
        nodes are removed. The dict represents an initial set of tags
        defining the permissible forward nodes from a given node (i.e.,
        those nodes lying on a cycle free path).
    """
    # Initialize the cycle-free paths graph and the tag dictionary
    dic_PG = {0: _initialize_cfpg(pg, source, target)}
    round_counter = 1
    # Perform CFPG generation in successive rounds to ensure convergence
    while True:
        #print("Starting round %d" % round_counter)
        #print("Level 0: %d nodes, %d edges" % (len(dic_PG[0][0]),
                                               #len(dic_PG[0][0].edges())))
        for k in range(1, path_length):
            # Start by copying the information from the previous level
            H = dic_PG[k-1][0].copy()
            tags = deepcopy(dic_PG[k-1][1])
            # Check if we have already detected there are no cycle free paths.
            # If so just propagate this information.
            if not H:
                dic_PG[k] = dic_PG[k-1]
            else:
                # Identify the nodes at level k in G_(k-1)
                X = [v for v in H.nodes_iter() if v[0] == k]
                # We will track the (g_x, tags_x) pairs contributed by each x
                # through dic_X
                dic_X = {}
                for x in X:
                    tags_x = {}
                    g_x_f = _forward(x, H, path_length)
                    g_x_b = _backward(x, H)
                    g_x = nx.DiGraph()
                    g_x.add_edges_from(g_x_b.edges())
                    g_x.add_edges_from(g_x_f.edges())
                    # Get the nodes in the forward reach set representing cycles
                    # back through node x, (excluding x at level k)
                    nodes_to_prune = [v for v in g_x_f
                                      if v[1] == x[1] and v[0] != k]
                    # If there are no nodes to prune then just add the tag 'x'
                    # to all the nodes in g_x_f but not to x
                    g_x_prune = _prune(g_x, nodes_to_prune, source, target)
                    # If target or x gets pruned then x will contribute
                    # nothing to G_k
                    if (target not in g_x_prune) or (x not in g_x_prune):
                        pass
                    nodes_to_tag = [v for v in g_x_prune.nodes()
                                    if v[0] >= k]
                    # Otherwise add the tag x to the nodes in the strict
                    # future of x. update dic_X
                    for v in g_x_prune.nodes_iter():
                        if v[0] >= k:
                            D = tags[v]
                            D.append(x)
                            tags_x[v] = D
                        else:
                            tags_x[v] = tags[v]
                    dic_X[x] = (g_x_prune, tags_x)
                # We can now piece together the pairs in dic_X to obtain (G_k,
                # tags_k)
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
            #print("Level %d: %d nodes, %d edges" % (k, len(dic_PG[k][0]),
                                                    #len(dic_PG[k][0].edges())))
        if not dic_PG[len(dic_PG)-1][0] or \
           set(dic_PG[0][0].edges()) == set(dic_PG[len(dic_PG)-1][0].edges()):
            break
        else:
            dic_PG = {0: dic_PG[k]}
        round_counter += 1
    return dic_PG


def _initialize_cfpg(pg, source, target):
    """Initialize cycle free paths graph data structures.

    Parameters
    ----------
    pg : networkx.DiGraph()
        "Raw" (contains cycles) paths graph as created by
        :py:func:`indra.explanation.paths_graph.paths_graph`.
    source : tuple
        Source node, of the form (0, source_name).
    target : tuple
        Target node, of the form (target_depth, source_name).

    Returns
    -------
    tuple : (networkx.DiGraph(), dict)
        The initialized, not-yet cycle free paths graph consists of the
        paths graph remaining after cycles through the source or target
        nodes are removed. The dict represents an initial set of tags
        defining the permissible forward nodes from a given node (i.e.,
        those nodes lying on a cycle free path).
    """
    # Identify the initial set of nodes to be pruned. In this initial phase,
    # they are simply nodes whose names match the source or target.
    nodes_to_prune = set([v for v in pg.nodes_iter()
                          if (v != source) and (v != target) and \
                             ((v[1] == source[1]) or (v[1] == target[1]))])
    # Get the paths graph after initial source/target cycle pruning
    pg_0 = _prune(pg, nodes_to_prune, source, target)
    # Initialize an empty list of tags for each node
    tags = dict([(node, []) for node in pg_0.nodes_iter()])
    # Add source tag to all nodes except source itself
    _add_tag(tags, source, [v for v in pg_0.nodes()])
    return (pg_0, tags)


def _prune(pg, nodes_to_prune, source, target):
    """Iteratively prunes nodes from a copy of the paths graph.

    We prune the graph *pg* iteratively by the following procedure:
      1. Remove the nodes given by *nodes_to_prune* from the graph.
      2. Identify nodes (other than the source node) that now have no
         incoming edges.
      3. Identify nodes (other than the target node) that now have no outgoing
         edges.
      4. Set *nodes_to_prune* to the nodes identified in steps 2 and 3.
      5. Repeat from 1 until there are no more nodes to prune.

    Parameters
    ----------
    pg : networkx.DiGraph
        Paths graph to prune.
    nodes_to_prune : list
        Nodes to prune from paths graph.
    source : tuple
        Source node, of the form (0, source_name).
    target : tuple
        Target node, of the form (target_depth, source_name).

    Returns
    -------
    networkx.DiGraph()
        Pruned paths graph.
    """
    # First check if we are pruning any nodes to prevent unnecessary copying
    # of the paths graph
    if not nodes_to_prune:
        return pg
    # Make a copy of the graph
    pg_pruned = pg.copy()
    # Perform iterative pruning
    while nodes_to_prune:
        # Remove the nodes in our pruning list
        pg_pruned.remove_nodes_from(nodes_to_prune)
        # Make a list of nodes whose in or out degree is now 0 (making
        # sure to exclude the source and target, whose depths are at 0 and
        # path_length, respectively)
        no_in_edges = [node for node, in_deg in pg_pruned.in_degree_iter()
                        if in_deg == 0 and node != source]
        no_out_edges = [node for node, out_deg in pg_pruned.out_degree_iter()
                        if out_deg == 0 and node != target]
        nodes_to_prune = set(no_in_edges + no_out_edges)
    return pg_pruned


# Function for updating node tags in place
def _add_tag(tag_dict, tag_node, nodes_to_tag):
    for v in nodes_to_tag:
        tag_dict[v].append(tag_node)


def _forward(v, H, length):
    """Compute the subgraph of H defined by the paths forward from node v.

    Parameters
    ----------
    v : tuple(int, str)
        The node to get the _forward subgraph for.
    H : networkx.DiGraph()
        For a given path length n, H defines the graph G_i at the i-th stage
        for 1 <= i <= n.
    Returns
    -------
    networkx.DiGraph()
        Subgraph reachable by forward paths from v in H.
    """
    j = v[0]
    L = {}
    L[j] = [v]
    h = nx.DiGraph()
    for k in range(j+1, length+1):
        for v in L[k-1]:
            h.add_edges_from(H.out_edges(v))
        L[k] = [w for w in h if w[0] == k]
    return h


def _backward(v, H):
    """Compute the subgraph of H defined by the paths backward from node v.

    Parameters
    ----------
    v : tuple(int, str)
        The node to get the _backward subgraph for.
    H : networkx.DiGraph()
        For a given path length n, H defines the graph G_i at the i-th stage
        for 1 <= i <= n.

    Returns
    -------
    networkx.DiGraph()
        Subgraph reachable by backward paths from v in H.
    """
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

def PG_cf(src, tgt, g, t):
    """

    Implements the major step (the outer loop) for constructing G_cf. We do so
    by computing dic_CF, a dictionary based version of G_cf.  dic_CF[i] will be
    a quadruple of the form (V_i, next_i, pred_i, t_i).

    V_i will be the set of nodes at level i.

    A node -after dic_CF[i] has been computed- will be of the form (i, n, c) where i is the level, n is the name and c is 
    the copy number of the node (i,n) in G_0. In other words, each node in G_0 will be split into 
    one or more copies to implement our memoryless sampling strategy.

    next_i is the successor relation

    pred_i[v] is the set of predecessors of of v in V_i. The construction proceeds from 
    the target to source. At stage i of the construction we convert nodes of the form (i, n) into 
    nodes of the form (i,n,c). For any such new node pred_i[v] will be nodes of the form 
    (i-1,n) at level i-1.

    t_i[v] will be the new tags of the node v. They will be pairs of the form (j,n).
    In other words their type will be the same as of T_0.
    (Note: In T_0,  I assign nodes of G_0 as tags rather than their names. This turns out 
    to be convenient for the construction of G_cf)

    Once the construction of PG_cf is complete we will no onger 
    require pred_i and t_i.
  

"""

    We first hardwire dic_CF[tgt[0]]. We just need to create one instance of tgt.
    """
    
    next_tgt = {}
    pred_tgt = {}
    tgt_0 = (tgt[0], tgt[1], 0)
    next_tgt[tgt_0] = []
    pred_tgt[tgt_0] = g.predecessors(tgt)
    t_cf_tgt = {}
    t_cf_tgt[tgt_0] = t[tgt]
    
    dic_CF = {tgt[0]: ([tgt_0], next_tgt, pred_tgt, t_cf_tgt)}
    
    """here is the main construction
    """
    empty = ([],{}, {},{})
    for i in reversed(range(1, tgt[0])):
        V_i_1, next_i_1, pred_i_1, t_cf_i_1 = dic_CF[i+1]
        
        if V_i_1 == []:
            for j in range(0, tgt[0]+1):
                dic_CF[j] = empty
            return dic_CF
        V_current = []
        for v in V_i_1:
            V_current.extend(pred_i_1[v])
            V_current = list(set(V_current))
        
        """
        Thus V_current is the set of nodes (which will be pairs)  at level i to be processed.
        the converted  nodes(which will be triples)  will be binned into V_i
        """
        if V_current == []:
            for j in range(0, tgt[0]+1):
                dic_CF[j] = empty
            return dic_CF
            
        else:
            V_i = []
            next_i = {}
            pred_i = {}
            t_cf_i = {}
            """ Now comes the heart of the construction. We take a node x in V_current
                 and split it into -in general- multiple copies to ensure that 
                 if (u,v) is an edge in G_cf then the set of tags of u is included in the
                 set of tags of v"""
            for x in V_current:
                """X_up is the set of nodes at the level above (remember tgt is above src!).
                These nodes will be triples.
                
                X_down is the nodes at the level below. They will be pairs
                """
                X_up = [w for w in V_i_1 if x in pred_i_1[w]]
                X_down = g.predecessors(x)
                if X_up == []:
                    continue
                """The actual splitting up and how to connect the resulting copies of x to its neighbors 
                above and below is carried out by the split_graph function. We explain this function later
                """
                
                V_x, next_x, pred_x, t_cf_x = split_graph(src, tgt, x,  X_up, X_down, t_cf_i_1, t, g)
                
                """We now extend V_i, next_i, pred_i and t_i in the obvious way.
                """
                
                V_i.extend(V_x)
                for u in V_x:
                    next_i[u] = next_x[u]
                    pred_i[u] = pred_x[u]
                    t_cf_i[u] = t_cf_x[u]
            dic_CF[i] = (V_i, next_i, pred_i, t_cf_i)               
    """Finally we hardwire dic_CF[0]
    """
         
    V_1 = dic_CF[1][0]
    src_0 = (src[0], src[1], 0)
    tgt_0 = (tgt[0], tgt[1], 0)
    V_0 = [src_0]
    next_src = {}    
    next_src[src_0] = V_1
    pred_src = {}
    pred_src[src_0] = []
    t_cf_src = {}
    t_cf_src[src_0] = t[src] 
    
    dic_CF[0] = (V_0 , next_src, pred_src, t_cf_src)        
    
    return dic_CF

def split_graph(src, tgt, x,  X_up, X_down, t_cf, t, g):
    V_x = []
    next_x = {}
    pred_x = {}
    t_x = {}
    S_up = {}
    
    
    """
    For each member of X_up we create one copy of x in terms of the tags it is going to get.
    This is the set X_wx below for each w in X_up. However X_wx is the set of nodes (in G_0) from which we can reach x
    without encountering x[1] AND without encountering w[1]. As result some nodes in X_wx may be isolated.
    Hence we prune them away. 
    """
    for w in X_up:
        X_wx =  set(t_cf[w]) & set(t[x])
        N_wx = list(X_wx)
        g_wx = g.subgraph(N_wx)
        nodes_prune = [v for v in g_wx if (v != x and g_wx.successors(v) == []) or (v!= src and g_wx.predecessors(v) == [])]
        g_wx_pruned = _prune(g_wx, nodes_prune, src, x)
        if x in g_wx_pruned and src in g_wx_pruned:
            s = frozenset(g_wx_pruned.nodes())
            S_up[w] = s
     
    
    """ Thus S_up[w] are the elements of S that point to w
    """
    S = []               
    for m in range(len(X_up)):
        w = X_up[m]
        if w in S_up.keys():
            S.append(S_up[w])
    S = list(set(S))
        
    """
each  element of S will be a pruned list of tags. We will creat one copy x_r
of x  for each element r in S. We assign r to be the set of tags of x_r.
next of x_r is assemebled using S_up and pred is defined in the expected way using X_down.
    """             
        
    c = 0
    for r in S:
        x_c = (x[0], x[1], c)
        V_x.append(x_c)
        next_x[x_c] = [w for w in X_up if w in S_up.keys() and r == S_up[w]]
        pred_x[x_c] = [u for u in X_down if u in r]
        t_x[x_c] = r
        c += 1

    return (V_x, next_x, pred_x, t_x)

def dic_to_graph(dic):
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


def dc_prev(src_0,tgt_0, dic):
    dic_prev = {}
    for k in range(len(dic.keys())):
        V_k, next_k,pred_k,t_k = dic[k]
        prev_k = {}
        if k == 0:
            prev_k[src_0] = []
            dic_prev[k] = (V_k, next_k, prev_k)
        else:
            V_k_1, next_k_1 ,pred_k_1, t_k_1 = dic[k-1]
            for v in V_k:
                prev_k[v] = [u for u in V_k_1 if v in next_k_1[u]]
            dic_prev[k] = (V_k, next_k, prev_k)
    return dic_prev
def PG_cf_pruned(src_0, tgt_0, dic):
    l = len(dic.keys())
    current = dic.copy()
    while True:
        I = isolated(src_0, tgt_0, current)
        if I == []:
            return current
        new = {}
        for i in range(l):
            V_i, next_i, prev_i = current[i]
            V_i_pruned = [v for v in V_i if v not in I]
            next_i_pruned = {}
            prev_i_pruned = {}
            
            for v in V_i_pruned:
                next_i_pruned[v]= [u for u in next_i[v] if u not in I]
                prev_i_pruned[v] = [u for u in prev_i[v] if u not in I]
                
            new[i] = (V_i_pruned, next_i_pruned, prev_i_pruned)
        current = new


def isolated(src_0, tgt_0, dic):
    l = len(dic.keys())
    isolated = []
    for i in range(1, l):
        V_i, next_i, prev_i = dic[i]
        V_i_1, next_i_1, prev_i_1 = dic[i-1]
        for v in V_i:
            prev_v = prev_i[v]
            if (prev_v == [] and v != src_0)  or (v != tgt_0 and  next_i[v] == []):
                isolated.append(v)
    isolated = list(set(isolated))
    return isolated

                
        
        
        
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

def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')

if __name__ == '__main__':
    """
    We use 25 randomly generated graphs for testing the algorithm
    """
    dic_rg_long = pickle.load(open("save.dic_rg_long", "rb"))
    for i in range(25):
        
        G_i, source, target = dic_rg_long[i] 
        #draw(G_i, 'G_1.pdf')
        
        print "graph#, no of nodes, no of edges", (i, len(G_i), len(G_i.edges()))
        
        (f_level, b_level)  =  paths_graph.get_reachable_sets(G_i, source, target,
                                              max_depth=14, signed=False)
        
        for length in range(5,10):
            print "paths length", length
            print "_______"
            
            P = list(nx.all_simple_paths(G_i, source, target, length+1))
            P_correct = [tuple(p) for p in P if len(p) == length+1]
            """
            For validation, we compute explicitly the set of paths in the original graph 
            of a fixed length
            """
            pg_raw = paths_graph.paths_graph(G_i, source, target, length, f_level,
                                     b_level, signed=False, target_polarity=0)
            src = (0, source)
            tgt = (length, target)
        
            dic_PG = cycle_free_paths_graph(pg_raw, src, tgt, length)
              
            G_0, T_0 = dic_PG[len(dic_PG)-1]
            """
            The above is the out of the iterative method
            """
            
            T_0[tgt].append(tgt)
            dic_cf = PG_cf(src,tgt,G_0,T_0)
            empty = ([],{}, {},{})
            if dic_cf[0] == empty:
                continue
            """dic_cf is the dictionary representation of the cycle free representation
            """
            src_0 = (src[0], src[1], 0)
            tgt_0 = (tgt[0], tgt[1], 0)
            G_cf = dic_to_graph(dic_cf)
            dic_cf_prev = dc_prev(src_0,tgt_0, dic_cf)
            
            """We throw away the pred_i and t_i dictionaries from dic_cf. We add the dictionary prev_i
            which yields for each v at levl i
            its set of predecessors in dic_cf
            """
            dic_cf_pruned = PG_cf_pruned(src_0, tgt_0, dic_cf_prev)
            
            """In general there will be isolated nodes in the cycle free representation.
            They don't seem to cause any trouble in terms of sampling for paths.
            This seems to suggest they are always isolated from below. No matter.
            We prune them away.
        
            We then convert the modified and reduced representation into a digraph.
            """
            G_cf_pruned = dic_to_graph(dic_cf_pruned)
            
            """We verify the three required properties.
            Recall:
            CF1: Every source-to-target path in G_cf is cycle free.
            CF2: Every cycle free path in the original graph appears as a source-to-target path 
                 in G_cf.
            CF3: There is a 1-1 correspondence between the paths in G_cf and the paths in the 
                 original graph. This means there is no redundancy in the representation.
                 For every path in the original graph there is a unique path in G_cf that corresponds to it.
            
            To do so, we first compute the set of source-to-target paths
            (the nodes will be triples) in G_cf_pruned
            """
            P_cf_pruned = list(nx.all_simple_paths(G_cf_pruned, src_0, tgt_0))
            
            """Next we extract the actual paths by projecting down to second component.
            """
            P_cf_pruned_names = name_paths(P_cf_pruned)
            
            """We first verify CF1.
            """
            
            for p in P_cf_pruned_names:
                if len(p) != len(list(set(p))):
                    print "cycle!"
                    print p
                    break 
            
            """Next we verify CF2. We will in fact check if the set of paths in
            P_cf_pruned_names is exactly the set of paths in the original graph.
            """
            if set(P_correct) != set(P_cf_pruned_names):
                print "Missing"
                print "graph, length", (i, length)
                print "______________"
            
            """Finally we verify CF3
            """
            
            if len(P_cf_pruned) != len(list(set(P_cf_pruned_names))):
                print "redundant representation!"
                print "graph, length", (i, length)
                print "____________"
                    
                    
                            
                        
                
            
            
           
