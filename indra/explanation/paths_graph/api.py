
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
u
