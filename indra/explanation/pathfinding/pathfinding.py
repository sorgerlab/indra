__all__ = ['shortest_simple_paths', 'bfs_search', 'find_sources',
           'get_path_iter', 'bfs_search_multiple_nodes']
import sys
import logging
from collections import deque
from copy import deepcopy

import networkx as nx
import networkx.algorithms.simple_paths as simple_paths

from heapq import heappush, heappop
from itertools import count

from numpy import log as ln

from .util import get_sorted_neighbors


logger = logging.getLogger(__name__)


# Copy from networkx.algorithms.simple_paths
# Added ignore_nodes and ignore_edges arguments
def shortest_simple_paths(G, source, target, weight=None, ignore_nodes=None,
                          ignore_edges=None, hashes=None, ref_counts_function=None, strict_mesh_id_filtering=False,
                          const_c=1, const_tk=10):
    """Generate all simple paths in the graph G from source to target,
       starting from shortest ones.

    A simple path is a path with no repeated nodes.

    If a weighted shortest path search is to be used, no negative weights
    are allowed.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    target : node
       Ending node for path

    weight : string
        Name of the edge attribute to be used as a weight. If None all
        edges are considered to have unit weight. Default value None.

    ignore_nodes : container of nodes
       nodes to ignore, optional

    ignore_edges : container of edges
       edges to ignore, optional

    hashes : list
        hashes specifying (if not empty) allowed edges
    
    ref_count_function : function
        function counting references and PMIDs of an edge from its
        statement hashes

    const_c : int
            Constant used in MeSH IDs-based weight calculation

    const_tk : int
        Constant used in MeSH IDs-based weight calculation

    Returns
    -------
    path_generator: generator
       A generator that produces lists of simple paths, in order from
       shortest to longest.

    Raises
    ------
    NetworkXNoPath
       If no path exists between source and target.

    NetworkXError
       If source or target nodes are not in the input graph.

    NetworkXNotImplemented
       If the input graph is a Multi[Di]Graph.

    Examples
    --------

    >>> G = nx.cycle_graph(7)
    >>> paths = list(nx.shortest_simple_paths(G, 0, 3))
    >>> print(paths)
    [[0, 1, 2, 3], [0, 6, 5, 4, 3]]

    You can use this function to efficiently compute the k shortest/best
    paths between two nodes.

    >>> from itertools import islice
    >>> def k_shortest_paths(G, source, target, k, weight=None):
    ...     return list(islice(nx.shortest_simple_paths(G, source, target,
    ...         weight=weight), k))
    >>> for path in k_shortest_paths(G, 0, 3, 2):
    ...     print(path)
    [0, 1, 2, 3]
    [0, 6, 5, 4, 3]

    Notes
    -----
    This procedure is based on algorithm by Jin Y. Yen [1]_.  Finding
    the first $K$ paths requires $O(KN^3)$ operations.

    See Also
    --------
    all_shortest_paths
    shortest_path
    all_simple_paths

    References
    ----------
    .. [1] Jin Y. Yen, "Finding the K Shortest Loopless Paths in a
       Network", Management Science, Vol. 17, No. 11, Theory Series
       (Jul., 1971), pp. 712-716.

    """
    if source not in G:
        s = source[0] if isinstance(source, tuple) else source
        raise nx.NodeNotFound('source node %s not in graph' % s)

    if target not in G:
        t = target[0] if isinstance(target, tuple) else target
        raise nx.NodeNotFound('target node %s not in graph' % t)

    if hashes:
        if strict_mesh_id_filtering:
            length_func = len
            shortest_path_func = _bidirectional_shortest_path
        else:
            weight = 'context_weight'
            def length_func(path):
                return sum(G.adj[u][v][weight] for (u, v) in zip(path, path[1:]))
            def shortest_path_func(G, source, target, weight, ignore_nodes, ignore_edges,
                                   force_edges):
                return simple_paths._bidirectional_dijkstra(G, source, target, weight,
                                                            ignore_nodes, ignore_edges)
            collect_weights = lambda path : [G[u][v]['context_weight'] for u, v in zip(path[:-1], path[1:])]
    else:
        if strict_mesh_id_filtering:
            return []
        if weight is None:
            length_func = len
            shortest_path_func = _bidirectional_shortest_path
            collect_weights = lambda path : ['N/A'] * (len(path) - 1)
        else:
            def length_func(path):
                return sum(G.adj[u][v][weight] for (u, v) in zip(path, path[1:]))
            def shortest_path_func(G, source, target, weight, ignore_nodes, ignore_edges,
                                   force_edges):
                return simple_paths._bidirectional_dijkstra(G, source, target, weight,
                                                            ignore_nodes, ignore_edges)
            collect_weights = lambda path : [G[u][v][weight] for u, v in zip(path[:-1], path[1:])]

    allowed_edges = []
    if hashes:
        for u, v, data, in G.edges(data=True):
            ref_counts, total = ref_counts_function([d['stmt_hash'] for d in data['statements']])
            if not ref_counts:
                ref_counts = 1e-15
            data['context_weight'] = -const_c * ln(ref_counts / (total + const_tk))

    culled_ignored_nodes = set() if ignore_nodes is None else set(ignore_nodes)
    culled_ignored_edges = set() if ignore_edges is None else set(ignore_edges)
    listA = list()
    listB = simple_paths.PathBuffer()
    prev_path = None
    while True:
        cur_ignore_nodes = culled_ignored_nodes.copy()
        cur_ignore_edges = culled_ignored_edges.copy()
        if not prev_path:
            length, path = shortest_path_func(G, source, target, weight=weight,
                                              ignore_nodes=cur_ignore_nodes,
                                              ignore_edges=cur_ignore_edges,
                                              force_edges=allowed_edges)
            listB.push(length, path)
        else:
            for i in range(1, len(prev_path)):
                root = prev_path[:i]
                root_length = length_func(root)
                for path in listA:
                    if path[:i] == root:
                        cur_ignore_edges.add((path[i - 1], path[i]))
                try:
                    length, spur = shortest_path_func(
                        G, root[-1], target, ignore_nodes=cur_ignore_nodes,
                        ignore_edges=cur_ignore_edges, weight=weight, 
                        force_edges=allowed_edges)
                    path = root[:-1] + spur
                    listB.push(root_length + length, path)
                except nx.NetworkXNoPath:
                    pass
                cur_ignore_nodes.add(root[-1])
        if listB:
            path = listB.pop()
            rcvd_ignore_values = yield path, collect_weights(path)
            if rcvd_ignore_values is not None:
                culled_ignored_nodes = culled_ignored_nodes.union(
                    rcvd_ignore_values[0])
                culled_ignored_edges = culled_ignored_edges.union(
                    rcvd_ignore_values[1])
            listA.append(path)
            prev_path = path
        else:
            break


# Implementation inspired by networkx's
# networkx.algorithms.traversal.breadth_first_search::generic_bfs_edges
def bfs_search(g, source_node, reverse=False, depth_limit=2, path_limit=None,
               max_per_node=5, node_filter=None, node_blacklist=None,
               terminal_ns=None, sign=None, max_memory=int(2**29), hashes=None, 
               strict_mesh_id_filtering=False, **kwargs):
    """Do breadth first search from a given node and yield paths

    Parameters
    ----------
    g : nx.Digraph
        An nx.DiGraph to search in. Can also be a signed node graph. It is
        required that node data contains 'ns' (namespace) and edge data
        contains 'belief'.
    source_node : node
        Node in the graph to start from.
    reverse : bool
        If True go upstream from source, otherwise go downstream. Default:
        False.
    depth_limit : int
        Stop when all paths with this many edges have been found. Default: 2.
    path_limit : int
        The maximum number of paths to return. Default: no limit.
    max_per_node : int
        The maximum number of paths to yield per parent node. If 1 is
        chosen, the search only goes down to the leaf node of its first
        encountered branch. Default: 5
    node_filter : list[str]
        The allowed namespaces (node attribute 'ns') for the nodes in the
        path
    node_blacklist : set[node]
        A set of nodes to ignore. Default: None.
    terminal_ns : list[str]
        Force a path to terminate when any of the namespaces in this list
        are encountered and only yield paths that terminate at these
        namepsaces
    sign : int
        If set, defines the search to be a signed search. Default: None.
    max_memory : int
        The maximum memory usage in bytes allowed for the variables queue
        and visited. Default: 1073741824 bytes (== 1 GiB).
    hashes : list
        List of hashes used (if not empty) to select edges for path finding

    Yields
    ------
    path : tuple(node)
        Paths in the bfs search starting from `source`.
    """
    int_plus = 0
    int_minus = 1

    allowed_edges = []
    if strict_mesh_id_filtering:
        if hashes:
            edge_by_hash = g.graph['edge_by_hash']
            for h in hashes:
                if h in edge_by_hash:
                    allowed_edges.append(edge_by_hash[h])
        else:
            return []

    queue = deque([(source_node,)])
    visited = ({source_node}).union(node_blacklist) \
        if node_blacklist else {source_node}
    yielded_paths = 0
    while queue:
        cur_path = queue.popleft()
        last_node = cur_path[-1]

        # if last node is in terminal_ns, continue to next path
        if terminal_ns and g.nodes[last_node]['ns'].lower() in terminal_ns \
                and source_node != last_node:
            continue

        sorted_neighbors = get_sorted_neighbors(G=g, node=last_node,
                                                reverse=reverse, 
                                                force_edges=allowed_edges)
        yielded_neighbors = 0
        # for neighb in neighbors:
        for neighb in sorted_neighbors:
            neig_name = neighb[0] if isinstance(neighb, tuple) else neighb

            # Check cycles
            if sign is not None:
                # Avoid signed paths ending up on the opposite sign of the
                # same node
                if (neig_name, int_minus) in cur_path or \
                        (neig_name, int_plus) in cur_path:
                    continue
            elif neighb in visited:
                continue

            # Check namespace
            if node_filter and len(node_filter) > 0:
                if g.nodes[neighb]['ns'].lower() not in node_filter:
                    continue

            # Add to visited nodes and create new path
            visited.add(neighb)
            new_path = cur_path + (neighb,)

            # Check yield and break conditions
            if len(new_path) > depth_limit + 1:
                continue
            elif not terminal_ns:
                # Yield newest path and recieve new ignore values

                # Signed search yield
                if sign is not None:
                    if reverse:
                        # Upstream signed search should not end in negative
                        # node
                        if new_path[-1][1] == int_minus:
                            ign_vals = None
                            pass
                        else:
                            ign_vals = yield new_path, ['N/A'] * (len(new_path) - 1)
                            yielded_paths += 1
                            yielded_neighbors += 1

                    else:
                        # Downstream signed search has to end on node with
                        # requested sign
                        if new_path[-1][1] != sign:
                            ign_vals = None
                            pass
                        else:
                            ign_vals = yield new_path, ['N/A'] * (len(new_path) - 1)
                            yielded_paths += 1
                            yielded_neighbors += 1

                # Unsigned search
                else:
                    ign_vals = yield new_path, ['N/A'] * (len(new_path) - 1)
                    yielded_paths += 1
                    yielded_neighbors += 1

            # terminal_ns is not None: only yield if last node is in
            # teminal_ns
            else:
                # If terminal_ns
                if g.nodes[neighb]['ns'].lower() in terminal_ns:
                    # If signed, reverse, negative start node OR
                    #    signed, not reverse, wrong sign:
                    # don't yield this path
                    if sign is not None and \
                            reverse and \
                            new_path[-1][1] == int_minus \
                            or \
                            sign is not None and \
                            not reverse and \
                            new_path[-1][1] != sign:
                        ign_vals = None
                        pass
                    else:
                        ign_vals = yield new_path, ['N/A'] * (len(new_path) - 1)
                        yielded_paths += 1
                        yielded_neighbors += 1
                else:
                    ign_vals = None
                    pass

            # If new ignore nodes are recieved, update set
            if ign_vals is not None:
                ign_nodes, ign_edges = ign_vals
                visited.update(ign_nodes)

            # Check max paths reached, no need to add to queue
            if path_limit and yielded_paths >= path_limit:
                break

            # Append yielded path
            queue.append(new_path)

            # Check for memory
            if sys.getsizeof(queue) + sys.getsizeof(visited) > max_memory:
                logger.warning('Memory overflow reached: %d' %
                               (sys.getsizeof(queue) + sys.getsizeof(visited)))
                raise StopIteration('Reached maxmimum allowed memory usage')

            # Check if we've visited enough neighbors
            # Todo: add all neighbors to 'visited' and add all skipped
            #  paths to queue? Currently only yielded paths are
            #  investigated deeper
            if max_per_node and yielded_neighbors >= max_per_node:
                break

        # Check path limit again to catch the inner break for path_limit
        if path_limit and yielded_paths >= path_limit:
            break


def bfs_search_multiple_nodes(g, source_nodes, path_limit=None, **kwargs):
    """Do breadth first search from each of given nodes and yield paths
    until path limit is met.

    Parameters
    ----------
    g : nx.Digraph
        An nx.DiGraph to search in. Can also be a signed node graph. It is
        required that node data contains 'ns' (namespace) and edge data
        contains 'belief'.
    source_nodes : list[node]
        List of nodes in the graph to start from.
    path_limit : int
        The maximum number of paths to return. Default: no limit.
    **kwargs : keyword arguments
        Any kwargs to pass to bfs_search.

    Yields
    ------
    path : tuple(node)
        Paths in the bfs search starting from `source`.
    """
    yielded_paths = 0
    for n in source_nodes:
        paths = bfs_search(g, n, path_limit=path_limit, **kwargs)
        for p in paths:
            yield p
            yielded_paths += 1
            if path_limit and yielded_paths >= path_limit:
                break
        if path_limit and yielded_paths >= path_limit:
            break


def get_path_iter(graph, source, target, path_length, loop, dummy_target):
    """Return a generator of paths with path_length cutoff from source to
    target."""
    path_iter = nx.all_simple_paths(graph, source, target, path_length)
    try:
        for p in path_iter:
            path = deepcopy(p)
            # Remove common target from a path.
            if dummy_target:
                path.remove(target)
            if loop:
                path.append(path[0])
            # A path should contain at least one edge
            if len(path) < 2:
                continue
            yield path
    except nx.NetworkXNoPath:
        pass


def find_sources(graph, target, sources):
    """Get the set of source nodes with paths to the target.

    Given a common target and  a list of sources (or None if test statement
    subject is None), perform a breadth-first search upstream from the
    target to determine whether there are any sources that have paths to
    the target. For efficiency, does not return the full path,
    but identifies the upstream sources and the length of the path.

    Parameters
    ----------
    graph : nx.DiGraph
        A DiGraph with signed nodes to find paths in.
    target : node
        The signed node (usually common target node) in the graph to start
        looking upstream for matching sources.
    sources : list[node]
        Signed nodes corresponding to the subject or upstream influence
        being checked.

    Returns
    -------
    generator of (source, path_length)
        Yields tuples of source node and path length (int). If there are no
        paths to any of the given source nodes, the generator is empty.
    """
    # First, create a list of visited nodes
    # Adapted from
    # networkx.algorithms.traversal.breadth_first_search.bfs_edges
    visited = set([target])
    # Generate list of predecessor nodes with a sign updated according to
    # the sign of the target node

    # The queue holds tuples of "parents" (in this case downstream nodes)
    # and their "children" (in this case their upstream influencers)
    queue = deque([(target, graph.predecessors(target), 0)])
    while queue:
        parent, children, path_length = queue[0]
        try:
            # Get the next child in the list
            child = next(children)
            # Is this child one of the source nodes we're looking for? If
            # so, yield it along with path length.
            # Also make sure that found source is positive
            if (sources is None or child in sources) and child[1] == 0:
                logger.debug("Found path to %s from %s with length %d"
                             % (target, child, path_length+1))
                yield (child, path_length+1)
            # Check this child against the visited list. If we haven't
            # visited it already (accounting for the path to the node),
            # then add it to the queue.
            if child not in visited:
                visited.add(child)
                queue.append(
                    (child, graph.predecessors(child), path_length + 1))
        # Once we've finished iterating over the children of the current
        # node, pop the node off and go to the next one in the queue
        except StopIteration:
            queue.popleft()
    # There was no path; this will produce an empty generator
    return

def _bidirectional_shortest_path(G, source, target,
                                 ignore_nodes=None,
                                 ignore_edges=None,
                                 weight=None,
                                 hashes=None,
                                 force_edges=None):
    """Returns the shortest path between source and target ignoring
       nodes and edges in the containers ignore_nodes and ignore_edges.

    This is a custom modification of the standard bidirectional shortest
    path implementation at networkx.algorithms.unweighted

    Parameters
    ----------
    G : NetworkX graph

    source : node
       starting node for path

    target : node
       ending node for path

    ignore_nodes : container of nodes
       nodes to ignore, optional

    ignore_edges : container of edges
       edges to ignore, optional

    weight : None
       This function accepts a weight argument for convenience of
       shortest_simple_paths function. It will be ignored.

    force_edges : list
        list specifying (if not empty) allowed edges

    Returns
    -------
    path: list
       List of nodes in a path from source to target.

    Raises
    ------
    NetworkXNoPath
       If no path exists between source and target.

    See Also
    --------
    shortest_path

    """
    # call helper to do the real work
    results = _bidirectional_pred_succ(G, source, target, ignore_nodes,
                                       ignore_edges, force_edges=force_edges)
    pred, succ, w = results

    # build path from pred+w+succ
    path = []
    # from w to target
    while w is not None:
        path.append(w)
        w = succ[w]
    # from source to w
    w = pred[path[0]]
    while w is not None:
        path.insert(0, w)
        w = pred[w]

    return len(path), path


def _bidirectional_pred_succ(G, source, target, ignore_nodes=None,
                             ignore_edges=None, force_edges=None):
    """Bidirectional shortest path helper.
       Returns (pred,succ,w) where
       pred is a dictionary of predecessors from w to the source, and
       succ is a dictionary of successors from w to the target.
    """
    # does BFS from both source and target and meets in the middle
    if ignore_nodes and (source in ignore_nodes or target in ignore_nodes):
        raise nx.NetworkXNoPath("No path between %s and %s."
                                % (source, target))
    if target == source:
        return ({target: None}, {source: None}, source)

    # handle either directed or undirected
    if G.is_directed():
        Gpred = G.predecessors
        Gsucc = G.successors
    else:
        Gpred = G.neighbors
        Gsucc = G.neighbors

    # support optional nodes filter
    if ignore_nodes:
        def filter_iter(nodes):
            def iterate(v):
                for w in nodes(v):
                    if w not in ignore_nodes:
                        yield w
            return iterate

        Gpred = filter_iter(Gpred)
        Gsucc = filter_iter(Gsucc)

    # support optional edges filter
    if ignore_edges or force_edges:
        if G.is_directed():
            def filter_pred_iter(pred_iter):
                def iterate(v):
                    if force_edges:
                        for w in pred_iter(v):
                            if (w, v) not in ignore_edges and (w, v)\
                                in force_edges:
                                    yield w
                    else:
                        for w in pred_iter(v):
                            if (w, v) not in ignore_edges:
                                yield w
                return iterate

            def filter_succ_iter(succ_iter):
                def iterate(v):
                    if force_edges:
                        for w in succ_iter(v):
                            if (v, w) not in ignore_edges and (v, w)\
                                in force_edges:
                                    yield w
                    else:
                        for w in succ_iter(v):
                            if (v, w) not in ignore_edges:
                                yield w
                return iterate

            Gpred = filter_pred_iter(Gpred)
            Gsucc = filter_succ_iter(Gsucc)

        else:
            def filter_iter(nodes):
                def iterate(v):
                    if force_edges:
                        for w in nodes(v):
                            if (v, w) not in ignore_edges \
                                and (w, v) not in ignore_edges \
                                    and (v, w) in force_edges and (w, v)\
                                        in force_edges:
                                            yield w
                    else:
                        for w in nodes(v):
                            if (v, w) not in ignore_edges \
                                and (w, v) not in ignore_edges:
                                    yield w
                return iterate

            Gpred = filter_iter(Gpred)
            Gsucc = filter_iter(Gsucc)

    # predecesssor and successors in search
    pred = {source: None}
    succ = {target: None}

    # initialize fringes, start with forward
    forward_fringe = [source]
    reverse_fringe = [target]

    while forward_fringe and reverse_fringe:
        if len(forward_fringe) <= len(reverse_fringe):
            this_level = forward_fringe
            forward_fringe = []
            for v in this_level:
                for w in Gsucc(v):
                    if w not in pred:
                        forward_fringe.append(w)
                        pred[w] = v
                    if w in succ:
                        # found path
                        return pred, succ, w
        else:
            this_level = reverse_fringe
            reverse_fringe = []
            for v in this_level:
                for w in Gpred(v):
                    if w not in succ:
                        succ[w] = v
                        reverse_fringe.append(w)
                    if w in pred:
                        # found path
                        return pred, succ, w

    raise nx.NetworkXNoPath("No path between %s and %s." % (source, target))


def open_dijkstra_search(g, start, reverse=False, depth_limit=2,
                         path_limit=None, hashes=None, terminal_ns=None,
                         weight=None, ref_counts_function=None, const_c=1,
                         const_tk=10):
    """Do Dijkstra search from a given node and yield paths

    Parameters
    ----------
    g : nx.Digraph
        An nx.DiGraph to search in.
    start : node
        Node in the graph to start from.
    reverse : bool
        If True go upstream from source, otherwise go downstream. Default:
        False.
    depth_limit : int
        Stop when all paths with this many edges have been found. Default: 2.
    path_limit : int
        The maximum number of paths to return. Default: no limit.
    terminal_ns : list[str]
        Force a path to terminate when any of the namespaces in this list
        are encountered and only yield paths that terminate at these
        namepsaces
    hashes : list
        List of hashes used to set edge weights
    weight : str
        Name of edge's attribute used as its weight
    ref_count_function : function
        function counting references and PMIDs of an edge from its
        statement hashes
    const_c : int
        Constant used in MeSH IDs-based weight calculation
    const_tk : int
        Constant used in MeSH IDs-based weight calculation
    

    Yields
    ------
    path : tuple(node)
        Paths in the bfs search starting from `source`.
    """
    allowed_edges = []
    if hashes:
        edge_by_hash = g.graph['edge_by_hash']
        for h in hashes:
            if h in edge_by_hash:
                allowed_edges.append(edge_by_hash[h])
        for u, v, data in g.edges(data=True):
            ref_counts, total = ref_counts_function([d['stmt_hash']
                for d in data['statements']])
            if not ref_counts:
                ref_counts = 1e-15
            data['context_weight'] = -const_c * ln(ref_counts / (total + const_tk))
        collect_weights = lambda path : [g[u][v]['context_weight'] for u, v in zip(path[:-1], path[1:])]

    if reverse:
        g = g.reverse(copy=False)

    if terminal_ns: # If not set, terminal_ns will be an empty list []
        def proper_path(path):
            if g.nodes[path[-1]]['ns'].lower() not in terminal_ns:
                return False
            for n in path[:-1]:
                if g.nodes[n]['ns'].lower() in terminal_ns:
                    return False
            return True
    else:
        proper_path = lambda x: True
    paths = list(nx.single_source_dijkstra_path(g, start,
        weight=('context_weight' if hashes else weight)).values())[1:]
    if path_limit is not None:
        for p in paths:
            path_limit -= 1
            if proper_path(p):
                yield p, collect_weights(p)
            if not path_limit:
                break
    else:
        for p in paths:
            if proper_path(p):
                yield p, collect_weights(p)
