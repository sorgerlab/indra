import os
import numpy as np
import networkx as nx
from indra import has_config

class PathsTree(object):
    """Build a tree representing a set of paths.

    Nodes in the tree are tuples representing the common prefix of all
    downstream paths. The head of the tree is an empty tuple, `()`. Each leaf
    of the tree represents a complete path.

    Parameters
    ----------
    paths : iterable of tuples
        Each element of the iterable is a tuple representing a sequence of
        nodes that constitutes a path.
    source_graph : networkx.DiGraph (optional)
        Source graph used to generate the set of paths and containing edge
        weights keyed by 'weight' in the edge metadata. If provided, allows
        weighted sampling over paths in the PathsTree. If not provided, all
        edges are considered to have locally equal weights of 1.

    Attributes
    ----------
    graph : networkx.DiGraph
        A directed graph representing the set of paths as a tree.
    """
    def __init__(self, paths, source_graph=None):
        self.graph = nx.DiGraph()
        if paths:
            edge_set = set()
            for path in paths:
                # Split path at all branch points
                for i in range(0, len(path)):
                    head = tuple(path[0:i])
                    tail = tuple(path[0:i+1])
                    edge_set.add((head, tail))
            # Get edge weights (have to do this as a separate step because the
            # weight dictionary is not hashable in the set
            edges_with_weights = []
            for head, tail in edge_set:
                if source_graph and len(head) > 0:
                    u = head[-1]
                    v = tail[-1]
                    weight = source_graph[u][v].get('weight', 1)
                else:
                    weight = 1
                edges_with_weights.append((head, tail, {'weight': weight}))
            self.graph.add_edges_from(edges_with_weights)

    def sample(self, num_samples=1000):
        """Sample a set of paths from the path tree.

        At each sampling step, the next node is chosen at random from the set
        of successors of the current node according to the edge weight in the
        'weight' entry of the edge data.

        Parameters
        ----------
        num_samples : int
            Number of paths to sample.
        """
        # Make sure we have a graph to sample from
        if not self.graph:
            return []
        # If so, do the sampling
        sampled_paths = []
        while len(sampled_paths) < num_samples:
            # The root of the tree should be the empty tuple
            node = tuple()
            while True:
                out_edges = self.graph.out_edges(node, data=True)
                # If there are no successors to the current node, then we've
                # hit a leaf of the tree and have found a path
                if not out_edges:
                    break
                # For determinism in testing, sort the out edges
                if has_config('TEST_FLAG'):
                    out_edges.sort()
                # The float is necessary here for Python 2 compatibility
                weights = [float(t[2]['weight']) for t in out_edges]
                # Normalize the weights to a proper probability distribution
                p = np.array(weights) / np.sum(weights)
                # Choose a successor at random based on the weights
                pred_idx = np.random.choice(range(len(out_edges)), p=p)
                node = out_edges[pred_idx][1]
            # Add the path (contained by the leaf node) to the list of sampled
            # paths
            sampled_paths.append(node)
        return sampled_paths



