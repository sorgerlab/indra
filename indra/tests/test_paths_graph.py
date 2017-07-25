import networkx
import pygraphviz
from indra.explanation import paths_graph

graph1 = networkx.DiGraph()
graph1.add_nodes_from(['A', 'B', 'C', 'D'])
graph1.add_edges_from([('A', 'B', {'polarity': 0}), ('B', 'D', {'polarity': 0}),
                      ('A', 'C', {'polarity': 0}), ('C', 'D', {'polarity': 0})])

graph2 = networkx.DiGraph()
graph2.add_nodes_from(['A', 'B', 'C', 'D'])
graph2.add_edges_from([('A', 'B', {'polarity': 0}), ('B', 'A', {'polarity': 0}),
                      ('A', 'C', {'polarity': 0}), ('C', 'D', {'polarity': 0})])

source = 'A'
target = 'D'
target_polarity = 0

def test_get_reachable_sets():
    f_level, b_level = paths_graph.get_reachable_sets(graph1, source, target)
    assert f_level == {0: {('A', 0)}, 1: {('B', 0), ('C', 0)}, 2: {('D', 0)}}
    assert b_level == {0: {('D', 0)}, 1: {('B', 0), ('C', 0)}, 2: {('A', 0)}}

if __name__ == '__main__':
    path_length = 1
    f_level, b_level = paths_graph.get_reachable_sets(graph1, source, target)
    pg = paths_graph.paths_graph(graph1, source, target, target_polarity,
                                 path_length, f_level, b_level)
    # No paths of length 1
    assert len(pg) == 0
    path_length = 2
    pg = paths_graph.paths_graph(graph1, source, target, target_polarity,
                                 path_length, f_level, b_level)
    paths = list(networkx.shortest_simple_paths(pg, (2, ('A', 0)), (0, ('D', 0))))
    assert len(paths) == 2
    assert [(2, ('A', 0)), (1, ('C', 0)), (0, ('D', 0))] in paths
    assert [(2, ('A', 0)), (1, ('B', 0)), (0, ('D', 0))] in paths
    #ag = networkx.nx_agraph.to_agraph(pg)
    #ag.draw('paths_graph_%d.pdf' % path_length, prog='dot')

