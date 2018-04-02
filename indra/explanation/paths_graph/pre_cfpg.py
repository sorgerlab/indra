import os
import logging
from copy import copy, deepcopy
import numpy as np
import networkx as nx
from .paths_graph import PathsGraph, PathSamplingException
from indra import has_config

logger = logging.getLogger('pre_cfpg')


class PreCFPG(PathsGraph):
    """Representation of a pre- cycle free paths graph with associated methods.

    The pre- cycle free paths graph consists of the paths graph remaining after
    cycles through the source or target nodes are removed. However, paths
    through the pre-CFPG node structure itself are not guaranteed to be cycle
    free; instead, cycle-free paths can be sampled by taking into account
    the tags associated with each node, representing the possible cycle-free
    histories of the node in terms of other upstream nodes.

    As with the "raw" paths graph (containing cycles), nodes in the pre-CFPG
    consist of tuples with two elements: (depth, name).

    Starting from the "raw" (i.e., containing cycles) paths graph, and
    given a target path length n, the algorithm iterates over each "level"
    in the graph 0 <= k <= n where level 0 consists only of the source node
    and level n consists only of the target.

    Each level k consists of a set of nodes, X; we examine each node x in X
    and identify the subset of nodes that are reachable in both the forward
    and backward directions from x. If any of the nodes in the forward
    reach subgraph contain x itself (but at a different depth), this
    represents a cyclic path back through x that is then pruned.

    Each node x therefore defines its own subgraph of cycle free paths,
    g_x.  After iterating over all x in X, we combine these subgraphs into
    the (in-progress) cycle free paths graph H_k. H_k therefore consists of
    the superset of nodes of all the subgraphs g_x for level k. When
    merging these subgraphs we prevent the re-introduction of cyclic paths
    by annotating each node in the graph with a list of "tags". The tags
    for any given node consist of a list of nodes lying at prior (upstream)
    levels. Therefore during sampling, transitions from an upstream node to
    a downstream node are only permissible if all nodes in the path up to a
    certain level are contained in the tag set of the downstream node.

    Parameters
    ----------
    pg : PathsGraph
        "Raw" (contains cycles) paths graph as created by
        :py:func:`indra.explanation.paths_graph.PathsGraph.from_graph`.
    graph : networkx.DiGraph
        The graph structure of the pre-CFPG.
    tags : dict
        A dictionary, keyed by node, with lists of other nodes representing
        the nodes lying upstream on cycle free paths. Node that each node
        also has itself as a tag.

    Attributes
    ----------
    source_node : tuple
        Node in the pre-CFPG graph representing the source: (0, source_name)
    target_node: tuple
        Node in the pre-CFPG graph representing the target:
        (path_length, target_name)
    """
    def __init__(self, pg, graph, tags):
        self.source_name = pg.source_name
        self.source_node = pg.source_node
        self.target_name = pg.target_name
        self.target_node = pg.target_node
        self.path_length = pg.path_length
        self.graph = graph
        self.tags = tags

    @classmethod
    def from_graph(klass, *args, **kwargs):
        """Compute a pre- cycle free paths graph from a graph.

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
        PreCFPG
            A instance of the PreCFPG the containing the pre- cycle free paths
            graph.
        """
        pg = PathsGraph.from_graph(*args, **kwargs)
        return PreCFPG.from_pg(pg)

    @classmethod
    def from_pg(klass, pg):
        """Compute a pre- cycle free paths graph from a PathsGraph.

        Parameters
        ----------
        pg : PathsGraph
            "Raw" (contains cycles) paths graph as created by
            :py:func:`indra.explanation.paths_graph.PathsGraph.from_graph`.

        Returns
        -------
        PreCFPG
            A instance of the PreCFPG the containing the pre- cycle free paths
            graph.
        """
        # Initialize the cycle-free paths graph and the tag dictionary
        source_node = pg.source_node
        target_node = pg.target_node
        dic_PG = {0: _initialize_pre_cfpg(pg)}
        round_counter = 1
        # Perform CFPG generation in successive rounds to ensure convergence
        logger.info("Creating pre-CFPG from PG")
        while True:
            logger.info("Starting round %d" % round_counter)
            for k in range(1, pg.path_length+1):
                logger.info("Iterating over level %d" % k)
                # Start by copying the information from the previous level
                H = dic_PG[k-1][0]
                tags = dic_PG[k-1][1]
                # Check if we have already detected there are no cycle free
                # paths, which would be indicated by an empty graph at the
                # previous level.  If so just propagate this information.
                if not H:
                    dic_PG[k] = dic_PG[k-1]
                else:
                    # Identify the nodes at level k in G_(k-1)
                    logger.info("Finding nodes at level %d" % k)
                    X = [v for v in H.nodes_iter() if v[0] == k]
                    # We will track the (g_x, tags_x) pairs contributed by each
                    # x through dic_X
                    dic_X = {}
                    logger.info("Iterating over nodes in level %d" % k)
                    for x in X:
                        tags_x = {}
                        g_x_f = _forward(x, H, pg.path_length)
                        g_x_b = _backward(x, H)
                        g_x = nx.DiGraph()
                        g_x.add_edges_from(g_x_b.edges(data=True))
                        g_x.add_edges_from(g_x_f.edges(data=True))
                        # Get the nodes in the forward reach set representing
                        # cycles back through node x, (excluding x at level k)
                        nodes_to_prune = [v for v in g_x_f
                                          if v[1] == x[1] and v[0] != k]

                        # If there are no nodes to prune then just add the tag
                        # 'x' to all the nodes in g_x_f but not to x
                        g_x_prune = prune(g_x, nodes_to_prune, source_node,
                                          target_node)
                        nodes_to_tag = [v for v in g_x_prune.nodes()
                                        if v[0] >= k]
                        # Otherwise add the tag x to the nodes in the strict
                        # future of x and update dic_X
                        for v in g_x_prune.nodes_iter():
                            if v[0] >= k:
                                D = tags[v]
                                D.append(x)
                                tags_x[v] = D
                            else:
                                tags_x[v] = tags[v]
                        dic_X[x] = (g_x_prune, tags_x)
                    # We can now piece together the pairs in dic_X to obtain
                    # (G_k, tags_k)
                    H_k = nx.DiGraph()
                    tags_k = {}
                    for x in X:
                        h_x = dic_X[x][0]
                        H_k.add_edges_from(h_x.edges(data=True))
                    # For every node in the combined subgraphs
                    for v in H_k.nodes_iter():
                        # Create a set of tags...
                        t = []
                        # ...by iterating over every node at this level
                        for x in X:
                            # ...and checking to see if the node v in the
                            # subgraph is in the history of a particular node x
                            # at this level
                            if v in dic_X[x][0]:
                                # ...if so, add v to the list of tags for x
                                tags_x = dic_X[x][1]
                                t.extend(tags_x[v])
                        t = list(set(t))
                        tags_k[v] = t
                    dic_PG[k] = (H_k, tags_k)
            if not dic_PG[len(dic_PG)-1][0] or \
               set(dic_PG[0][0].edges()) == \
                                set(dic_PG[len(dic_PG)-1][0].edges()):
                break
            else:
                dic_PG = {0: dic_PG[k]}
            round_counter += 1
        pre_cfpg, tags = dic_PG[pg.path_length]
        # Return the fully processed cfpg as an instance of the PreCFPG class
        return klass(pg, pre_cfpg, tags)

    def enumerate_paths(self):
        raise NotImplementedError()

    def count_paths(self):
        raise NotImplementedError()

    def set_uniform_path_distribution():
        raise NotImplementedError()

    def _successor(self, path, u):
        """Randomly choose a successor node of u given the current path."""
        out_edges = []
        for edge in self.graph.out_edges(u, data=True):
            if set(path) <= set(self.tags[edge[1]]):
                out_edges.append(edge)
        # If there are no admissible successors, raise a PathSamplingException
        if not out_edges:
            raise PathSamplingException("No cycle-free successors")
        # For determinism in testing
        if has_config('TEST_FLAG'):
            out_edges.sort()
        weights = [t[2]['weight'] for t in out_edges]
        # Normalize the weights to a proper probability distribution
        p = np.array(weights) / np.sum(weights)
        pred_idx = np.random.choice(len(out_edges), p=p)
        return out_edges[pred_idx][1]


def _initialize_pre_cfpg(pg):
    """Initialize pre- cycle free paths graph data structures.

    Parameters
    ----------
    pg : PathsGraph
        "Raw" (contains cycles) paths graph as created by
        :py:func:`indra.explanation.paths_graph.paths_graph`.
    source_node : tuple
        Source node, of the form (0, source_name).
    target_node : tuple
        Target node, of the form (target_depth, source_name).

    Returns
    -------
    tuple : (networkx.DiGraph(), dict)
    """
    # Identify the initial set of nodes to be pruned. In this initial phase,
    # they are simply nodes whose names match the source or target.
    nodes_to_prune = set([v for v in pg.graph.nodes_iter()
                        if (v != pg.source_node) and (v != pg.target_node) and \
                             ((v[1] == pg.source_node[1]) or \
                              (v[1] == pg.target_node[1]))])
    # Get the paths graph after initial source_node/target_node cycle pruning
    pre_cfpg_0 = prune(pg.graph, nodes_to_prune, pg.source_node, pg.target_node)
    # Initialize an empty list of tags for each node
    tags = dict([(node, []) for node in pre_cfpg_0.nodes_iter()])
    # Add source_node tag to all nodes
    _add_tag(tags, pg.source_node, [v for v in pre_cfpg_0.nodes()])
    return (pre_cfpg_0, tags)


def _add_tag(tag_dict, tag_node, nodes_to_tag):
    for v in nodes_to_tag:
        tag_dict[v].append(tag_node)


def prune(pg, nodes_to_prune, source, target):
    """Iteratively prunes nodes from (a copy of) a paths graph or CFPG.

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
        for v in L[k - 1]:
            h.add_edges_from(H.out_edges(v, data=True))
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
            h.add_edges_from(H.in_edges(v, data=True))
        L[k] = [w for w in h if w[0] == k]
    return h


