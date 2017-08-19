import networkx as nx
from indra.explanation import paths_graph as pg
from indra.explanation import cycle_free_paths as cfp

g1_uns = nx.DiGraph()
g1_uns.add_nodes_from(['A', 'B', 'C', 'D'])
g1_uns.add_edges_from((('A', 'B'), ('B', 'C'), ('C', 'D')))

g2_uns = nx.DiGraph()
g2_uns.add_nodes_from(['A', 'B', 'C', 'D'])
g2_uns.add_edges_from((('A', 'B'), ('B', 'A'), ('B', 'C'),
                      ('C', 'D'), ('A', 'D')))
source = 'A'
target = 'D'
length = 3

def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')

def test_prune():
    g = nx.DiGraph()
    g.add_edges_from((('S', 'A'), ('S', 'B'), ('A', 'S'), ('B', 'C'),
                      ('C', 'D'), ('D', 'T'), ('B', 'T')))
    length = 4
    (f_level, b_level) = pg.get_reachable_sets(g, 'S', 'T', max_depth=length)
    pg_raw = pg.paths_graph(g, 'S', 'T', length, f_level, b_level)
    pg_raw_edges = pg_raw.edges()
    nodes_to_prune = [(2, 'S')]
    # Prune the graph
    pg_pruned = cfp.prune(pg_raw, nodes_to_prune, length)
    # Make sure we didn't change the original graphs or node lists
    assert nodes_to_prune == [(2, 'S')]
    assert pg_raw.edges() == pg_raw_edges
    # The correctly pruned structure
    assert set(pg_pruned.edges()) == \
           set([((0, 'S'), (1, 'B')), ((1, 'B'), (2, 'C')),
                ((2, 'C'), (3, 'D')), ((3, 'D'), (4, 'T'))])

def test_pg_0():
    """
    # We first run the pg_0 calculation on a simple graph with no cycles
    # involving the source or target
    (f_level, b_level) = pg.get_reachable_sets(g1_uns, source, target,
                                               max_depth=length)
    pg_raw = pg.paths_graph(g1_uns, source, target, length, f_level, b_level)
    (pg_0, tags) = cfp.PG_0((0, source), (length, target), pg_raw)
    # Because no nodes are pruned, the initialized "cycle free"
    # paths graph will be the same as the path graph we started
    # with
    assert pg_0 == pg_raw
    assert tags == {(0, 'A'): [], (1, 'B'): ['A'], (2, 'C'): ['A'],
                    (3, 'D'): ['A']}
    """
    # The next graph contains a cycle passing through the source node, A
    (f_level, b_level) = pg.get_reachable_sets(g2_uns, source, target,
                                               max_depth=length)
    pg_raw = pg.paths_graph(g2_uns, source, target, length, f_level, b_level)
    (pg_0, tags) = cfp.PG_0((0, source), (length, target), pg_raw)


if __name__ == '__main__':
    test_prune()
    #test_pg_0()
