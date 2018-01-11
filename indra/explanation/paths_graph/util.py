

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

