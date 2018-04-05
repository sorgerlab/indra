import os
import itertools
import numpy as np
import networkx as nx
from indra import logging
from collections import defaultdict
from indra import has_config

logger = logging.getLogger('paths_graph')


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
    signed : boolean
        Whether the graph is signed. If True, sign information should be encoded
        in the 'sign' field of the edge data, with 0 indicating a positive edge
        and 1 indicating a negative edge.

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
        respectively. Note that if the source or target are not reachable by
        any path within the given maximum depth, both dicts are empty.
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
            return ({}, {})
    return (f_level, b_level)


class PathsGraph(object):
    def __init__(self, source_name, target_name, graph, path_length, signed,
                 target_polarity):
        self.source_name = source_name
        self.target_name = target_name
        self.signed = signed
        self.target_polarity = target_polarity
        if signed:
            self.source_name = (source_name, 0)
            self.source_node = (0, self.source_name)
            self.target_name = (target_name, target_polarity)
            self.target_node = (path_length, self.target_name)
        else:
            self.source_node = (0, source_name)
            self.target_node = (path_length, target_name)
        self.graph = graph
        self.path_length = path_length

    @classmethod
    def from_graph(klass, g, source, target, length, fwd_reachset=None,
                   back_reachset=None, signed=False, target_polarity=0):
        """Create a graph where all nodes lie on a path of the given length.

        Nodes in the graph account for the cumulative polarity from the source
        to the target, so that any path found from the source to the target
        will have the overall polarity as specified by the target_polarity
        argument.

        The function uses the forward and backward reach sets provided as
        arguments to efficiently identify the subset of nodes that are
        reachable from both the forward and backward directions. Pairs of nodes
        at neighboring levels are then checked against the original graph to
        determine which nodes are connected by edges with the appropriate
        direction and polarity. These nodes are then used to create a new
        graph, the "paths graph," which consists solely of these nodes and
        edges. This graph represents the superset of all possible paths from
        source to target of a given legnth and target polarity. Specific paths
        can then be obtained by sampling.

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
        fwd_reachset : Optional[dict]
            Dictionary of sets representing the forward reachset computed over
            the original graph g up to a maximum depth greater than the
            requested path length.  If not provided, the forward reach set is
            calculated up to the requested path length up to the requested path
            length by calling paths_graph.get_reachable_sets.
        back_reachset : Optional[dict]
            Dictionary of sets representing the backward reachset computed over
            the original graph g up to a maximum depth greater than the
            requested path length.  If not provided, the backward reach set is
            calculated up to the requested path length up to the requested path
            length by calling paths_graph.get_reachable_sets.
        signed : bool
            Specifies whether the underlying graph and the corresponding
            f_level and b_level reachable sets have signed edges.  If True,
            sign information should be encoded in the 'sign' field of the edge
            data, with 0 indicating a positive edge and 1 indicating a negative
            edge.
        target_polarity : 0 or 1
            Specifies the polarity of the target node: 0 indicates
            positive/activation, 1 indicates negative/inhibition.

        Returns
        -------
        PathsGraph
            Instance of PathsGraph class representing paths from source to
            target with a given length and overall polarity.
        """
        # If the reachable sets aren't provided by the user, compute them here
        # with a maximum depth given by the target path length.
        if fwd_reachset is None or back_reachset is None:
            (fwd_reachset, back_reachset) = get_reachable_sets(g, source,
                                                 target, max_depth=length,
                                                 signed=signed)

        # If either fwd_reachset or back_reachset is an empty dict (as they
        # would be if the nodes were unreachable from either directions) return
        # an empty paths graph
        if not fwd_reachset or not back_reachset:
            paths_graph = nx.DiGraph()
            return PathsGraph(source, target, paths_graph, length, signed,
                              target_polarity)
        # Otherwise, if the reachable sets are provided, use them after checking
        # if they have a depth at least equal to the given path length
        _check_reach_depth('forward', fwd_reachset, length)
        _check_reach_depth('backward', back_reachset, length)
        # Also, if the reachable sets do not have entries at the given length,
        # this means that either we are attempting to create a paths_graph for
        # a path longer than we generated reachable sets, or there is no path of
        # the given length (may depend on whether cycles were eliminated when
        # when generating the reachable sets).
        if not (length in fwd_reachset and length in back_reachset):
            paths_graph = nx.DiGraph()
            return klass(source, target, paths_graph, length, signed,
                         target_polarity)
        # By default, we set the "adjusted backward reach set", aka
        # back_reachset_adj, to be the same as the original back_reachset; this
        # is only overriden if we have a signed graph and a negative target
        # polarity
        back_reachset_adj = back_reachset
        # Signed graphs
        if signed:
            level = {0: set([(source, 0)]),
                     length: set([(target, target_polarity)])}
            # If the target polarity is even (positive/neutral), then the
            # cumulative polarities in the forward direction will match those
            # in the reverse direction; if the target polarity is odd, then the
            # polarities will be opposite at each matching node. Thus we check
            # the target polarity and flip the polarities of the backward reach
            # set if appropriate.
            if target_polarity == 1:
                back_reachset_adj = {}
                for i in range(0, len(back_reachset)):
                    polar_set = set()
                    for (u, w) in back_reachset[i]:
                        w_flipped = (w + 1) % 2
                        polar_set.add((u, w_flipped))
                    back_reachset_adj[i] = polar_set
        # Unsigned graphs
        else:
            level = {0: set([source]), length: set([target])}
        # Next we calculate the subset of nodes at each level that are reachable
        # from both the forward and backward directions. Because the polarities
        # have already been set appropriately, we can do this with a simple
        # set intersection.
        for i in range(1, length):
            f_reach_set = fwd_reachset[i]
            b_reach_set = back_reachset_adj[length - i]
            path_nodes = set(f_reach_set) & set(b_reach_set)
            level[i] = path_nodes
        # Next we explicitly enumerate the path graph nodes by tagging each node
        # with its level in the path
        pg_nodes = {}
        for i in range(0,length+1):
            pg_nodes[i] = list(itertools.product([i], level[i])) 
        # Finally we add edges between these nodes if they are found in the
        # original graph. Note that we have to check for an edge of the
        # appropriate polarity.
        pg_edges = []
        g_edges = []
        edge_weights = {}
        # Collect edge and edge weight info from the graph
        for u, v, data in g.edges_iter(data=True):
            if signed:
                edge_key = (u, v, data['sign'])
            else:
                edge_key = (u, v)
            g_edges.append(edge_key)
            edge_weights[edge_key] = float(data.get('weight', 1.0))
        for i in range(0, length):
            actual_edges = []
            logger.info("paths_graph: identifying edges at level %d" % i)
            if signed:
                # This set stores the information for performing the set
                # intersection with the edges in the source graph
                possible_edges = set()
                # This dict stores the information for the actual edge, with
                # weight, as we will need to add it to the PG
                edge_lookup = {}
                for edge in itertools.product(pg_nodes[i], pg_nodes[i+1]):
                    u_name, u_pol = edge[0][1]
                    v_name, v_pol = edge[1][1]
                    # If the polarity between neighboring nodes is the same,
                    # then we need a positive edge
                    required_sign = 0 if u_pol == v_pol else 1
                    edge_key = (u_name, v_name, required_sign)
                    possible_edges.add(edge_key)
                    edge_lookup[edge_key] = edge
                for edge_key in possible_edges.intersection(g_edges):
                    weighted_edge = edge_lookup[edge_key] + \
                                       ({'weight': edge_weights[edge_key]},)
                    actual_edges.append(weighted_edge)
            else:
                # Build a set representing possible edges between adjacent
                # levels
                possible_edges = set([(u[1], v[1]) for u, v in
                                itertools.product(pg_nodes[i], pg_nodes[i+1])])
                # Actual edges are the ones contained in the original graph; add
                # to list with prepended depths
                for u, v in possible_edges.intersection(g_edges):
                    actual_edges.append(((i, u), (i+1, v),
                                         {'weight': edge_weights[(u, v)]}))
            pg_edges.extend(actual_edges)
            logger.info("Done.")
        paths_graph = nx.DiGraph()
        paths_graph.add_edges_from(pg_edges)
        logger.info("Paths graph for length %d has %d nodes" %
                    (length, len(paths_graph)))
        return klass(source, target, paths_graph, length, signed,
                     target_polarity)

    def enumerate_paths(self, names_only=True):
        if not self.graph:
            return tuple([])
        paths = [tuple(path) for path in nx.all_simple_paths(self.graph,
                                    self.source_node, self.target_node)]
        if names_only:
            paths = self._name_paths(paths)
        return tuple(paths)

    def _get_path_counts(self):
        """Get a dictionary giving the number of paths through each node.

        The entry for the source node gives the total number of paths in the
        graph.
        """
        if not self.graph:
            return {}
        # Group nodes by level
        levels = defaultdict(list)
        for node in self.graph.nodes():
            levels[node[0]].append(node)
        # Initialize the path count
        path_counts = {}
        path_counts[self.target_node] = 1
        # Iterate over the levels going "backwards" from the target. This way
        # the path count at each node reflects the number of paths passing
        # through that node from source to target
        for i in reversed(range(0, self.path_length)):
            # Iterate over the nodes at this level
            for node in levels[i]:
                # The count for this node is the sum of the counts over all
                # of its successors
                path_counts[node] = \
                        sum([path_counts[succ]
                            for succ in self.graph.successors(node)])
        return path_counts

    def count_paths(self):
        """Count the total number of paths without enumerating them.

        Returns
        -------
        int
            The number of paths.
        """
        path_counts = self._get_path_counts()
        total_count = path_counts.get(self.source_node)
        if total_count is None:
            total_count = 0
        return total_count

    def set_uniform_path_distribution(self):
        """Adjusts edge weights to allow uniform sampling of paths.

        Note that calling this method will over-write any existing edge
        weights in the graph.
        """
        path_counts = self._get_path_counts()
        weight_dict = {}
        for u in self.graph.nodes():
            count_tuples = [(v, float(path_counts[v]))
                            for v in self.graph.successors(u)]
            if not count_tuples:
                continue
            v_list, counts = zip(*count_tuples)
            weights = np.array(counts) / np.sum(counts)
            for ix, v in enumerate(v_list):
                weight_dict[(u, v)] = weights[ix]
        nx.set_edge_attributes(self.graph, 'weight', weight_dict)

    @staticmethod
    def _name_paths(paths):
        return [tuple([node[1] for node in path]) for path in paths]

    def sample_paths(self, num_samples, names_only=True):
        """Sample paths of the given length between source and target.

        Parameters
        ----------
        num_samples : int
            The number of paths to sample.
        names_only : boolean
            Whether the paths should consist only of node names, or of node
            tuples (e.g., including depth and polarity). Default is True
            (only names).

        Returns
        -------
        list of tuples
            Each item in the list is a tuple of strings representing a path.
            Note that the paths may not be unique.
        """
        if not self.graph:
            return tuple([])
        paths = []
        while len(paths) < num_samples:
            try:
                path = self.sample_single_path(names_only=False)
                paths.append(path)
            except PathSamplingException:
                pass
        if names_only:
            paths = self._name_paths(paths)
        return tuple(paths)

    def sample_single_path(self, names_only=True):
        """Sample a path between source and target.

        Parameters
        ----------
        names_only : boolean
            Whether the paths should consist only of node names, or of node
            tuples (e.g., including depth and polarity). Default is True
            (only names).

        Returns
        -------
        tuple
            Tuple of nodes or node names representing a path.
        """
        # Sample a path from the paths graph.
        # If the path graph is empty, there are no paths
        if not self.graph:
            return tuple([])
        path = [self.source_node]
        current = self.source_node
        while current[1] != self.target_name:
            next = self._successor(path, current)
            path.append(next)
            current = next
        if names_only:
            path = tuple([node[1] for node in path])
        else:
            path = tuple(path)
        return path

    def _successor(self, path, node):
        out_edges = self.graph.out_edges(node, data=True)
        # For determinism in testing
        if has_config('TEST_FLAG'):
            out_edges.sort()
        weights = [t[2]['weight'] for t in out_edges]
        # Normalize the weights to a proper probability distribution
        p = np.array(weights) / np.sum(weights)
        pred_idx = np.random.choice(len(out_edges), p=p)
        return out_edges[pred_idx][1]


class CombinedPathsGraph(object):
    """Combine PathsGraphs for different lengths into a single super-PG.

    This is particularly useful for sampling paths where the sampled paths
    reflect the likelihood of drawing paths of particular lengths, given
    the weights on the edges.


    Parameters
    ----------
    pg_list : list of cfpg instances

    Attributes
    ----------
    source_name
    source_node
    target_name
    target_node
    graph
    """
    def __init__(self, pg_list):
        self.graph = nx.DiGraph()
        for pg in pg_list:
            self.graph.add_edges_from(pg.graph.edges(data=True))
        # Add info from the last PG in the list
        self.source_name = pg.source_name
        self.source_node = pg.source_node
        self.target_name = pg.target_name
        self.target_node = pg.target_node
        self.signed = pg.signed
        self.target_polarity = pg.target_polarity
        # Internally we create a PG wrapping the graph so as to re-use its
        # sampling method by composition
        self._pg = PathsGraph(pg.source_name, pg.target_name, self.graph,
                              None, pg.signed, pg.target_polarity)

    def sample_paths(self, num_samples):
        return self._pg.sample_paths(num_samples=num_samples)


def _check_reach_depth(dir_name, reachset, length):
    depth = max(reachset.keys())
    if depth < length:
        logger.warning("Insufficient depth: path length is %d "
                       "but %s reach set has maximum depth %d " %
                       (length, dir_name, depth))


class PathSamplingException(Exception):
    """Indicates a problem with sampling, e.g. a dead-end in the Pre-CFPG."""
    pass

