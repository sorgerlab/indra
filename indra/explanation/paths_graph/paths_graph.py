import os
import random
import itertools
import numpy as np
import networkx as nx
from indra import logging


logger = logging.getLogger('paths_graph')


def get_edges(sif_file):
    """Load edges from a SIF file with lines of the form 'u polarity v'"""
    edges = []
    with open(sif_file, 'rt') as f:
        for line in f.readlines():
            u, polarity, v = line.strip().split(' ')
            # Eliminate self-loops
            if u == v:
                pass
            else:
                edges.append((u, v, {'sign': int(polarity)}))
    g = nx.DiGraph()
    g.add_edges_from(edges)
    return g


def get_reachable_sets(g, source, target, max_depth=10, signed=False):
    """Get sets of nodes reachable from source and target at different depths.

    Parameters
    ----------
    g : nx.DiGraph
        The underlying graph used for computing reachable node sets.
    source : str
        Name of source node.
    target : target
        Name of target node.
    max_depth : int
        Maximum path length (depth) over which to compute reachable sets.

    Returns
    -------
    tuple
        The first item in the tuple represents the nodes reachable from the
        source in the forward direction; the second represents the nodes
        reachable from the target in the backward direction. Each reachable
        set takes the form of a dict with integers representing different
        depths as keys, and sets of nodes as values. The nodes themselves
        consist of tuples (u, w), where u is the name of the node and w is the
        cumulative polarity, forwards or backwards, from the source or target,
        respectively.
    """
    # Forward and backward level sets for signed and unsigned graphs
    if signed:
        source = (source, 0)
        target = (target, 0)
        f_level = {0: set([source])}
        b_level = {0: set([target])}
    else:
        f_level = {0: set([source])}
        b_level = {0: set([target])}
    # A bit of trickery to avoid a duplicated for loop--may be too much!
    directions = (
      ('forward', f_level, lambda u: [((u, v), v) for v in g.successors(u)]),
      ('backward', b_level, lambda v: [((u, v), u) for u in g.predecessors(v)]))
    # Utility function to make code below more compact
    def _add_signed_edge(reachable_set, node_polarity, edge_polarity):
        cum_polarity = (node_polarity + edge_polarity) % 2
        reachable_set.add((reachable_node, cum_polarity))
    # Iterate over levels
    for direction, level, edge_func in directions:
        visited = set([source]) if direction == 'forward' else set([target])
        for i in range(1, max_depth+1):
            reachable_set = set()
            # Signed graph
            if signed:
                for node, node_polarity in level[i-1]:
                    for (u, v), reachable_node in edge_func(node):
                        edge_dict = g.get_edge_data(u, v)
                        edge_polarity = edge_dict.get('sign')
                        # If this is a multidigraph, get_edge_data will return
                        # a dict keyed by integers
                        if edge_polarity is None:
                            for edge_key, edge_data in edge_dict.items():
                                _add_signed_edge(reachable_set, node_polarity,
                                                 edge_data['sign'])
                        else:
                            _add_signed_edge(reachable_set, node_polarity,
                                             edge_polarity)
            # Unsigned graph
            else:
                for node in level[i-1]:
                    for (u, v), reachable_node in edge_func(node):
                        reachable_set.add(reachable_node)
            visited = visited | reachable_set
            # If the reachable set is empty then we can stop
            if not reachable_set:
                break
            else:
                level[i] = reachable_set
        # If we're going forward we make sure we visited the target
        if (direction == 'forward' and target not in visited) or \
           (direction == 'backward' and source not in visited):
            return (None, None)
    return (f_level, b_level)


def from_graph(g, source, target, length, f_level, b_level,
               signed=False, target_polarity=0):
    """Generate a graph where all nodes lie on a path of the given length.

    Nodes in the graph account for the cumulative polarity from the source
    to the target, so that any path found from the source to the target will
    have the overall polarity as specified by the target_polarity argument.

    The function uses the forward and backward reach sets provided as arguments
    to efficiently identify the subset of nodes that are reachable from both
    the forward and backward directions. Pairs of nodes at neighboring levels
    are then checked against the original graph to determine which nodes are
    connected by edges with the appropriate direction and polarity. These
    nodes are then used to create a new graph, the "paths graph," which
    consists solely of these nodes and edges. This graph represents the
    superset of all possible paths from source to target of a given legnth
    and target polarity. Specific paths can then be obtained by sampling.

    Parameters
    ----------
    g : networkx.DiGraph
        The underlying graph on which paths will be generated.
    source : str
        Name of the source node.
    target : str
        Name of the target node.
    target_polarity : int
        Whether the desired path from source to target is positive (0)
        or negative (1).
    length : int
        Length of paths to compute.
    f_level : dict of sets
        Sets of nodes reachable from the source node at different depths.
        Generated by get_reachable_sets.
    b_level : dict of sets
        Sets of nodes reachable (backwards) from the target node at different
        depths. Generated by get_reachable_sets.
    signed : bool
        Specifies whether the underlying graph and the corresponding
        f_level and b_level reachable sets have signed edges.

    Returns
    -------
    nx.DiGraph
        Graph representing paths from source to target with a given length
        and overall polarity. If there are no paths of the specified length,
        returns an empty graph.
    """
    # If either f_level or b_level is None (as they would be if the nodes
    # were unreachable from either directions) return an empty graph
    if f_level is None or b_level is None:
        return nx.DiGraph()
    # Also, if the reachable sets do not have entries at the given length,
    # this means that either we are attempting to create a paths_graph for
    # a path longer than we generated reachable sets, or there is no path of
    # the given length (may depend on whether cycles were eliminated when
    # when generating the reachable sets).
    if not (length in f_level and length in b_level):
        return nx.DiGraph()
    # By default, we set the "adjusted backward reach set", aka b_level_adj,
    # to be the same as the original b_level; this is only overriden if we
    # have a signed graph and a negative target polarity
    b_level_adj = b_level
    # Signed graphs
    if signed:
        level = {0: set([(source, 0)]),
                 length: set([(target, target_polarity)])}
        # If the target polarity is even (positive/neutral), then the
        # cumulative polarities in the forward direction will match those in
        # the reverse direction; if the target polarity is odd, then the
        # polarities will be opposite at each matching node. Thus we check the
        # target polarity and flip the polarities of the backward reach set if
        # appropriate.
        if target_polarity == 1:
            b_level_adj = {}
            for i in range(0, len(b_level)):
                polar_set = set()
                for (u, w) in b_level[i]:
                    w_flipped = (w + 1) % 2
                    polar_set.add((u, w_flipped))
                b_level_adj[i] = polar_set
    # Unsigned graphs
    else:
        level = {0: set([source]), length: set([target])}
    # Next we calculate the subset of nodes at each level that are reachable
    # from both the forward and backward directions. Because the polarities
    # have already been set appropriately, we can do this with a simple
    # set intersection.
    for i in range(1, length):
        f_reach_set = f_level[i]
        b_reach_set = b_level_adj[length - i]
        path_nodes = set(f_reach_set) & set(b_reach_set)
        level[i] = path_nodes
    # Next we explicitly enumerate the path graph nodes by tagging each node
    # with its level in the path
    pg_nodes = {}
    for i in range(0,length+1):
        pg_nodes[i] = list(itertools.product([i], level[i])) 
    # Finally we add edges between these nodes if they are found in the original
    # graph. Note that we have to check for an edge of the appropriate polarity.
    pg_edges = set()
    if signed:
        g_edges_raw = g.edges(data=True)
        g_edges = [(u, v, data['sign']) for u, v, data in g_edges_raw]
    else:
        g_edges = [(u, v) for u, v in g.edges()]
    for i in range(0, length):
        actual_edges = set()
        logger.info("paths_graph: identifying edges at level %d" % i)
        if signed:
            possible_edges = set()
            edge_lookup = {}
            for edge in itertools.product(pg_nodes[i], pg_nodes[i+1]):
                u_name, u_pol = edge[0][1]
                v_name, v_pol = edge[1][1]
                # If the polarity between neighboring nodes is the same, then
                # we need a positive edge
                required_sign = 0 if u_pol == v_pol else 1
                edge_key = (u_name, v_name, required_sign)
                possible_edges.add(edge_key)
                edge_lookup[edge_key] = edge
            for edge_key in possible_edges.intersection(g_edges):
                actual_edges.add(edge_lookup[edge_key])
        else:
            # Build a set representing possible edges between adjacent levels
            possible_edges = set([(u[1], v[1]) for u, v in
                               itertools.product(pg_nodes[i], pg_nodes[i+1])])
            # Actual edges are the ones contained in the original graph; add
            # to list with prepended depths
            for u, v in possible_edges.intersection(g_edges):
                actual_edges.add(((i, u), (i+1, v)))
        pg_edges |= actual_edges
    path_graph = nx.DiGraph()
    path_graph.add_edges_from(pg_edges)
    return PathsGraph(source, target, path_graph, length, signed,
                      target_polarity)


class PathsGraph(object):
    def __init__(self, source_name, target_name, graph, path_length, signed,
                 target_polarity):
        self.source_name = source_name
        self.target_name = target_name
        self.signed = signed
        self.target_polarity = target_polarity
        if signed:
            self.source_node = (0, (source_name, 0))
            self.target_node = (path_length, (target_name, target_polarity))
        else:
            self.source_node = (0, source_name)
            self.target_node = (path_length, target_name)
        self.graph = graph
        self.path_length = path_length


    """
    def sample_single_path(self, weighted=False):
        # Sample a path from the paths graph.
        # If the path graph is empty, there are no paths
        if not self.graph:
            return tuple([])
        # Set the source node
        if self.signed:
            source_node = (0, (source, 0))
            target = (target, target_polarity)
        else:
            source_node = (0, source)
        assert source_node in pg
        # Repeat until we find a path without a cycle
        while True:
            path = [source_node[1]]
            current_node = source_node
            while True:
                if weighted:
                    out_edges = pg.out_edges(current_node, data=True)
                else:
                    out_edges = pg.out_edges(current_node)
                if 'TEST_FLAG' in os.environ:
                    out_edges.sort()
                if out_edges:
                    if weighted:
                        weights = [t[2]['weight'] for t in out_edges]
                        pred_idx = np.random.choice(range(len(out_edges)),
                                                    p=weights)
                        v = out_edges[pred_idx][1]
                    else:
                        v = np.random.choice(out_edges)[1]
                    # If we've already hit this node, it's a cycle; skip
                    if v[1] in path:
                        break
                    path.append(v[1])
                    if v[1] == target:
                        return tuple(path)
                current_node = v
    """

def combine_path_graphs(pg_dict):
    """Combine a dict of path graphs into a single super-pathgraph."""
    cpg = nx.DiGraph()
    for level, pg in pg_dict.items():
        # Start by adding
        for edge in pg:
            cpg.add_edges_from(pg.edges())
    return cpg


def sample_paths(g, source, target, max_depth=None, num_samples=1000,
                 eliminate_cycles=True, signed=False, target_polarity=0,
                 by_depth=False):
    # A helper function for running the sampling loop
    def _sample(pg, num):
        paths = []
        for i in range(num):
            path = sample_single_path(pg, source, target, signed=signed,
                                      target_polarity=target_polarity)
            if path:
                paths.append(path)
        return paths

    logger.info("Computing forward and backward reach sets...")
    # By default the max_depth is the number of nodes
    if max_depth is None:
        max_depth = len(g)
    f_level, b_level = get_reachable_sets(g, source, target, max_depth,
                                          signed=signed)
    # Compute path graphs over a range of path lengths
    pg_by_length = {}
    paths = []
    for path_length in range(1, max_depth+1):
        logger.info("Length %d: computing paths graph" % path_length)
        pg = paths_graph(g, source, target, path_length, f_level, b_level,
                         signed=signed, target_polarity=target_polarity)
        pg_by_length[path_length] = pg
        # If we're sampling by depth, do sampling here
        if pg and by_depth:
            logger.info("Length %d: Sampling %d paths" %
                        (path_length, num_samples))
            paths += _sample(pg, num_samples)
    # If we're sampling by depth, we've already collected all our paths
    if by_depth:
        return paths
    # Otherwise, run the sampling on the combined path graph
    else:
        # Combine the path graphs into one
        logger.info("Sampling %d paths from the combined path graph" %
                    num_samples)
        cpg = combine_path_graphs(pg_by_length)
        # If the combined path graph is empty, return an empty path list
        if not cpg:
            return []
        # Otherwise, sample from the combined path graph
        else:
            paths = _sample(pg, num_samples)
            return paths


def paths_to_graphset(paths_dict, pg_dict):
    from graphillion import GraphSet
    # Construct the universe
    edges = []
    nodes = set([])
    for path_length, pg in pg_dict.items():
        for e in pg.edges():
            # Putting in the depth seemed to cause problems
            #new_edge = tuple([str(v[1][0]) for v in e])
            new_edge = tuple([str(v[1]) for v in e])
            nodes |= set(new_edge)
            if new_edge in edges or tuple([new_edge[1], new_edge[0]]) in edges:
                continue
            else:
                edges.append(new_edge)
    # Set the graphset
    GraphSet.set_universe(edges)
    # Create a graphset from the paths in paths_dict
    all_paths = []
    for path_length, path_list in paths_dict.items():
        for path in path_list:
            new_path = []
            for node_ix in range(len(path) - 1):
                from_node = str(path[node_ix])
                to_node = str(path[node_ix + 1])
                assert from_node in nodes
                assert to_node in nodes
                new_path.append((from_node, to_node))
            all_paths.append(new_path)
    gs = GraphSet(all_paths)
    return gs


if __name__ == '__main__':
    g = get_edges('korkut_im.sif')
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'
    target_polarity = 0
    im_paths = set(sample_paths(g, source, target, signed=True,
                                target_polarity=0, by_depth=True, max_depth=8,
                                num_samples=100000))
    # 101 paths
    #gs = paths_to_graphset(paths_dict, pg_dict)

    g = get_edges('korkut_model_pysb_pysb.sif')
    source = 'BLK'
    target = 'EIF4EBP1'
    sif_paths = set(sample_paths(g, source, target, signed=False, by_depth=True,
                             max_depth=8, num_samples=10000000))
    #gs = paths_to_graphset(paths_dict, pg_dict)
    # 78,057 paths
