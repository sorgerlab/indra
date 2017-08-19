import networkx
from indra.explanation import paths_graph as pg
from indra.explanation import cycle_free_paths as cfp

g1_uns = networkx.DiGraph()
g1_uns.add_nodes_from(['A', 'B', 'C', 'D'])
g1_uns.add_edges_from((('A', 'B'), ('B', 'C'), ('C', 'D')))

g2_uns = networkx.DiGraph()
g2_uns.add_nodes_from(['A', 'B', 'C', 'D'])
g2_uns.add_edges_from((('A', 'B'), ('B', 'A'), ('B', 'C'),
                      ('C', 'D'), ('A', 'D')))
source = 'A'
target = 'D'
length = 3

def test_pg_0():
    (f_level, b_level) = pg.get_reachable_sets(
                                        g1_uns, source, target,
                                        max_depth=length)
    pg_raw = pg.paths_graph(g1_uns, source, target, length,
                            f_level, b_level)
    (pg_0, tags) = cfp.PG_0((0, source), (length, target), pg_raw)
    # Because no nodes are pruned, the initialized "cycle free"
    # paths graph will be the same as the path graph we started
    # with
    assert pg_0 == pg_raw
    assert tags == {(0, 'A'): [], (1, 'B'): ['A'],
                    (2, 'C'): ['A'], (3, 'D'): ['A']}

if __name__ == '__main__':
    test_pg_0()
