import random
import itertools
from indra import logging
import networkx as nx
import graphillion

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
                edges.append((u, v, {'polarity': int(polarity)}))
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
        Maximum depth over which to compute reachable sets.

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
    for direction, level, edge_func in directions:
        visited = set([source]) if direction == 'forward' else set([target])
        for i in range(1, max_depth):
            reachable_set = set()
            # Signed graph
            if signed:
                for node, node_polarity in level[i-1]:
                    for (u, v), reachable_node in edge_func(node):
                        edge_polarity = g.get_edge_data(u, v)['polarity']
                        cum_polarity = (node_polarity + edge_polarity) % 2
                        reachable_set.add((reachable_node, cum_polarity))
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


def paths_graph(g, source, target, length, f_level, b_level,
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
    for i in range(0, length):
        edges_at_this_level = set()
        for u, v in itertools.product(pg_nodes[i], pg_nodes[i+1]):
            # Signed graphs
            if signed:
                u_name, u_pol = u[1]
                v_name, v_pol = v[1]
                if (u_name, v_name) in g.edges():
                    edge_polarity = g.get_edge_data(u_name, v_name)['polarity']
                    # Look for an edge that flips or doesn't flip the polarity
                    # of the path depending on what we see in the cumulative
                    # polarities
                    if (u_pol == v_pol and edge_polarity == 0) or \
                       (u_pol != v_pol and edge_polarity == 1):
                        edges_at_this_level.add((u, v))
            # Unsigned graphs
            else:
                u_name = u[1]
                v_name = v[1]
                if (u_name, v_name) in g.edges():
                    edges_at_this_level.add((u, v))
        pg_edges |= edges_at_this_level
    path_graph = nx.DiGraph()
    path_graph.add_edges_from(pg_edges)
    return path_graph


def sample_single_path(g, source, target, length):
    """Sample a path from the paths graph."""
    path = []
    source_node = (length, (source, 0))
    path.append(source_node)
    assert source_node in g
    while True:
        current_node = path[-1]
        succs = g.successors(current_node)
        if succs:
            v = random.choice(succs)
            path.append(v)
        else:
            break
    # Remove the prepended node levels
    return tuple([u[1] for u in path])


def sample_paths(g, source, target, target_polarity, max_depth=10,
                 samples_per_length=100, draw_graphs=False):
    logger.info("Computing forward and backward reach sets...")
    f_level, b_level = get_reachable_sets(g, source, target, max_depth)
    # Compute path graph for a specific path length
    paths_by_length = {}
    pg_by_length = {}
    for path_length in range(1, max_depth + 1):
        logger.info("Length %d: computing paths graph" % path_length)
        pg = paths_graph(g, source, target, target_polarity, path_length,
                         f_level, b_level)

        # If there are no paths of this length, the graph will be empty
        if len(pg) == 0:
            logger.info("Length %d: no paths" % path_length)
            continue
        # Optionally draw the paths graph for paths of this length
        if draw_graphs:
            logger.info("Length %d: drawing Graphviz graph" % path_length)
            ag = nx.nx_agraph.to_agraph(pg)
            ag.draw('paths_graph_%d.pdf' % path_length, prog='dot')
        # Sample a given number of paths for this path length
        num_paths = 100
        logger.info("Length %d: Sampling %d paths" %
                    (path_length, samples_per_length))
        paths_for_length = []
        while len(paths_for_length) < num_paths:
            p = sample_single_path(pg, source, target, path_length)
            paths_for_length.append(p)
        # Remove redundant paths
        paths_by_length[path_length] = list(set(paths_for_length))
        pg_by_length[path_length] = pg
    return (paths_by_length, pg_by_length)


def paths_to_graphset(paths_dict, pg_dict):
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
    graphillion.GraphSet.set_universe(edges)
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
    gs = graphillion.GraphSet(all_paths)
    return gs


def combine_path_graphs(pg_dict):
    cpg = nx.DiGraph()
    for level, pg in pg_dict.items():
        # Start by adding
        for edge in pg:
            cpg.add_edges_from(pg.edges())
    return cpg

if __name__ == '__main__':
    g = get_edges('korkut_im.sif')
    source = 'BLK_phosphoY389_phosphorylation_PTK2_Y397'
    target = 'EIF4EBP1_T37_p_obs'
    target_polarity = 0
    logger.info("Testing the logger.")
    (paths_dict, pg_dict) = sample_paths(g, source, target, target_polarity,
                                         draw_graphs=True, max_depth=8)
    gs = paths_to_graphset(paths_dict, pg_dict)

