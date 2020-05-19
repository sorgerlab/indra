import networkx as nx
import networkx.algorithms.simple_paths as simple_paths


# Copy from networkx.algorithms.simple_paths
# Added ignore_nodes and ignore_edges arguments
def shortest_simple_paths(G, source, target, weight=None, ignore_nodes=None,
                          ignore_edges=None):
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

    if weight is None:
        length_func = len
        shortest_path_func = simple_paths._bidirectional_shortest_path
    else:
        def length_func(path):
            return sum(G.adj[u][v][weight] for (u, v) in zip(path, path[1:]))
        shortest_path_func = simple_paths._bidirectional_dijkstra

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
                                              ignore_edges=cur_ignore_edges)
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
                        ignore_edges=cur_ignore_edges, weight=weight)
                    path = root[:-1] + spur
                    listB.push(root_length + length, path)
                except nx.NetworkXNoPath:
                    pass
                cur_ignore_nodes.add(root[-1])
        if listB:
            path = listB.pop()
            rcvd_ignore_values = yield path
            if rcvd_ignore_values is not None:
                culled_ignored_nodes = culled_ignored_nodes.union(
                    rcvd_ignore_values[0])
                culled_ignored_edges = culled_ignored_edges.union(
                    rcvd_ignore_values[1])
            listA.append(path)
            prev_path = path
        else:
            break
