import logging
import itertools
import numpy as np
import networkx as nx
from .paths_graph import get_reachable_sets, PathsGraph
from .cfpg import CFPG

logger = logging.getLogger('paths_graph')


__all__ = ['load_signed_sif', 'sample_paths', 'enumerate_paths', 'count_paths']


def load_signed_sif(sif_file):
    """Load edges from a SIF file with lines of the form 'u polarity v'.

    Entries within each line can be separated by spaces and/or tabs. Polarity
    should be specified by 0 (for a positive/activating edge) or 1 (for a
    negative/inhibitory edge).

    Parameters
    ----------
    sif_file : str
        Path to the SIF file.

    Returns
    -------
    nx.DiGraph
        Graph with the sign information encoded in the 'sign' attribute of each
        edge.
    """
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


def sample_paths(g, source, target, max_depth=None, num_samples=1000,
                 cycle_free=True, signed=False, target_polarity=0):
    """Sample paths over a range of lengths from a graph.

    This high-level function provides explicit access to path sampling
    without the user need to explicit create PathsGraphs or CFPGs for
    different path lengths.

    Note that this function samples an equal number of paths from each depth;
    to sample paths where the sampling distribution reflects the probability
    of reach paths of different lengths, use an instance of `CombinedCFPG`.

    Parameters
    ----------
    g : networkx.DiGraph
        The underlying graph on which paths will be generated.
    source : str
        Name of the source node.
    target : str
        Name of the target node.
    max_depth : Optional[int]
        The maximum path length to consider. If not specified, the number of
        nodes in the graph is used as the default maximum depth.
    num_samples : int
        Number of path samples at each depth.
    cycle_free : bool
        If True, sample only cycle-free paths using CFPGs. Default is True.
    signed : bool
        Specifies whether the underlying graph and the corresponding
        f_level and b_level reachable sets have signed edges.  If True,
        sign information should be encoded in the 'sign' field of the edge
        data, with 0 indicating a positive edge and 1 indicating a negative
        edge.
    target_polarity : 0 or 1
        For a signed graph, specifies the polarity of the target node: 0
        indicates positive/activation, 1 indicates negative/inhibition.

    Returns
    -------
    list of paths
        Each path in the list contains a sequence of node names representing
        a path from source to target.
    """
    return _run_by_depth('sample_paths', [num_samples], g, source, target,
                         max_depth, cycle_free, signed, target_polarity)


def enumerate_paths(g, source, target, max_depth=None,
                    cycle_free=True, signed=False, target_polarity=0):
    """Enumerate paths over a range of lengths.

    Parameters
    ----------
    g : networkx.DiGraph
        The underlying graph on which paths will be generated.
    source : str
        Name of the source node.
    target : str
        Name of the target node.
    max_depth : Optional[int]
        The maximum path length to consider. If not specified, the number of
        nodes in the graph is used as the default maximum depth.
    cycle_free : bool
        If True, sample only cycle-free paths using CFPGs. Default is True.
    signed : bool
        Specifies whether the underlying graph and the corresponding
        f_level and b_level reachable sets have signed edges.  If True,
        sign information should be encoded in the 'sign' field of the edge
        data, with 0 indicating a positive edge and 1 indicating a negative
        edge.
    target_polarity : 0 or 1
        For a signed graph, specifies the polarity of the target node: 0
        indicates positive/activation, 1 indicates negative/inhibition.

    Returns
    -------
    list of paths
        Each path in the list contains a sequence of node names representing
        a path from source to target.
    """
    return _run_by_depth('enumerate_paths', [], g, source, target, max_depth,
                         cycle_free, signed, target_polarity)


def count_paths(g, source, target, max_depth=None,
                cycle_free=True, signed=False, target_polarity=0):
    """Count unique paths over a range of lengths without explicit enumeration.

    Parameters
    ----------
    g : networkx.DiGraph
        The underlying graph on which paths will be generated.
    source : str
        Name of the source node.
    target : str
        Name of the target node.
    max_depth : Optional[int]
        The maximum path length to consider. If not specified, the number of
        nodes in the graph is used as the default maximum depth.
    cycle_free : bool
        If True, sample only cycle-free paths using CFPGs. Default is True.
    signed : bool
        Specifies whether the underlying graph and the corresponding
        f_level and b_level reachable sets have signed edges.  If True,
        sign information should be encoded in the 'sign' field of the edge
        data, with 0 indicating a positive edge and 1 indicating a negative
        edge.
    target_polarity : 0 or 1
        For a signed graph, specifies the polarity of the target node: 0
        indicates positive/activation, 1 indicates negative/inhibition.

    Returns
    -------
    int
        Total number of paths up to the specified maximum depth.
    """
    return _run_by_depth('count_paths', [], g, source, target, max_depth,
                         cycle_free, signed, target_polarity)


def _run_by_depth(func_name, func_args, g, source, target, max_depth=None,
                  cycle_free=True, signed=False, target_polarity=0):
    """Run a function over paths graphs computed for different lengths."""
    if max_depth is None:
        max_depth = len(g)
    f_level, b_level = get_reachable_sets(g, source, target, max_depth,
                                          signed=signed)
    # Compute path graphs over a range of path lengths
    pg_by_length = {}
    if func_name == 'count_paths':
        results = 0
    else:
        results = []
    for path_length in range(1, max_depth+1):
        logger.info("Length %d: computing paths graph" % path_length)
        args = [g, source, target, path_length, f_level, b_level]
        kwargs = {'signed': signed, 'target_polarity': target_polarity}
        if cycle_free:
            pg = CFPG.from_graph(*args, **kwargs)
        else:
            pg = PathsGraph.from_graph(*args, **kwargs)
        pg_by_length[path_length] = pg
        # If we're sampling by depth, do sampling here
        if pg:
            func = getattr(pg, func_name)
            results += func(*func_args)
    return results

