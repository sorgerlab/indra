"""
Compute cycle-free paths graphs for use in sampling causal paths.

A *paths-graph,* say *G_n*, has the property that every node in *G_n* lies on a
path of length n from source to target. These paths may not be cycle-free. This
module transforms *G_n* into a new graph *G_cf* such that in *G_cf* every node
lies on a cycle-free path from the source to a target. Note that it will *not*
be the case, despite the misleading name *G_cf*, that *G_cf* contains *only*
cycle-free paths. However, we are able to "easily" sample cycle-free paths from
*G_cf* without any backtracking using metadata attached to each node in *G_cf*.

The algorithm is described further in :py:func:`cycle_free_paths_graph`, below.
"""

import numpy as np
import itertools
from copy import copy, deepcopy
import networkx as nx
from indra.explanation import paths_graph


def cycle_free_paths_graph(pg, source, target, path_length):
    """Compute a cycle free paths graph.

    Starting from the "raw" (i.e., containing cycles) paths graph, and given a
    target path length n, the algorithm iterates over each "level" in the
    graph 0 <= k <= n where level 0 consists only of the source node and level
    n consists only of the target.

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
        print("Starting round %d" % round_counter)
        print("Level 0: %d nodes, %d edges" % (len(dic_PG[0][0]),
                                               len(dic_PG[0][0].edges())))
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
                    g_x_f = _forward(x, H)
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
                    nodes_to_tag = [v for v in g_x_prune.nodes_iter()
                                    if v[0] > k]
                    # Otherwise add the tag x to the nodes in the strict
                    # future of x. update dic_X
                    for v in g_x_prune.nodes_iter():
                        if v[0] > k:
                            D = tags[v]
                            D.append(x[1])
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
            print("Level %d: %d nodes, %d edges" % (k, len(dic_PG[k][0]),
                                                    len(dic_PG[k][0].edges())))
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
    _add_tag(tags, source, [v for v in pg_0.nodes_iter() if v != source])
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
        tag_dict[v].append(tag_node[1])


def _forward(v, H):
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
    for k in range(j+1, 10):
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


def _cf_succ(H, t, path, v):
    """Randomly choose a successor node of v.

    Parameters
    ----------
    H : networkx.DiGraph()
        The cycle free paths graph.
    t : dict
        The tags dictionary.
    path : list
        The path so far (list of nodes).
    v : tuple
        The current node.

    Returns
    -------
    tuple
        Randomly chosen successor node on a non-cyclic path.
    """
    succ = []
    for u in H.successors(v):
        if set(path) <= set(t[u]):
            succ.append(u)
    # Note that the circuitous way of choosing from this list is the result of
    # the odd way numpy.random handles lists of lists (it excepts).
    idx_list = list(range(len(succ)))
    w_idx = np.random.choice(idx_list)
    w = succ[w_idx]
    return w


def sample_single_path(source, target, H, t):
    """Sample a single cycle-free path.

    The sampling procedure uses the tag sets to trace out cycle-free
    paths. If we have reached a node *v* via the path *p* then we can choose
    the successor *u* of *v* as the next node only if *p* appears in the tag
    set of u.

    Parameters
    ----------
    source : tuple
        Source node, of the form (0, source_name).
    target : tuple
        Target node, of the form (target_depth, source_name).
    H : networkx.DiGraph()
        The cycle free paths graph.
    t : dict
        The tags dictionary.

    Returns
    -------
    list of strings
        A randomly sampled, non-cyclic path. Nodes are represented as node
        names only, i.e., the depth prefixes are removed.
    """
    path = [source[1]]
    current = source
    while current != target:
        next = _cf_succ(H, t, path, current)
        """ a sanity check; since I have not stree-tested the code yet """
        assert next[1] not in path, "Error: found a cycle"
        path.append(next[1])
        current = next
    return tuple(path)


def sample_many_paths(source, target, H, t, n):
    """Sample many cycle-free paths.

    Parameters
    ----------
    source : tuple
        Source node, of the form (0, source_name).
    target : tuple
        Target node, of the form (target_depth, source_name).
    H : networkx.DiGraph()
        The cycle free paths graph.
    t : dict
        The tags dictionary.

    Returns
    -------
    list of lists
        Each item in the list is a list of strings representing a path. Note
        that the paths may not be unique.
    """
    # If the graph is empty, then there are no paths
    if not H:
        return []
    P = []
    for i in range(0, n):
        p = sample_single_path(source, target, H, t)
        P.append(p)
    return P


def sample_paths(g, source, target, max_depth, target_polarity=None,
                 num_paths=1000):
    print('Source: %s' % source)
    signed = True if target_polarity is not None else False

    (f_level, b_level)  = paths_graph.get_reachable_sets(
                                g, source, target, max_depth=10, signed=signed)
    all_paths = []
    for length in range(1, max_depth):
        print("Path length: %d" % length)
        pg_raw = paths_graph.paths_graph(g, source, target, length, f_level,
                             b_level, signed=signed,
                             target_polarity=target_polarity)
        # Append depths to our source and target nodes
        src = (0, (source, 0))
        tgt = (length, (target, target_polarity))
        dic_PG = cycle_free_paths_graph(pg_raw, src, tgt, length)
        G_cf, T = dic_PG[length-1]
        try:
            P = sample_many_paths(src, tgt, G_cf, T, num_paths)
            all_paths.extend(P)
        except IndexError:
            pass
    return all_paths


if __name__ == '__main__':
    G_0 = paths_graph.get_edges('korkut_im.sif')
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'

    (f_level, b_level)  =  paths_graph.get_reachable_sets(G_0, source, target,
                                              max_depth=10, signed=False)
    length = 8

    pg_raw = paths_graph.paths_graph(G_0, source, target, length, f_level,
                                     b_level, signed=False, target_polarity=0)

    # Append depths to our source and target nodes
    src = (0, source)
    tgt = (length, target)
    dic_PG = cycle_free_paths_graph(pg_raw, src, tgt, length)
    G_cf, T = dic_PG[7]
    P = sample_many_paths(src, tgt, G_cf, T, 1000)
    #print("--- %s seconds ---" % (time.time() - start_time))

