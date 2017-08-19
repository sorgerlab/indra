import time
import networkx as nx
from indra.explanation import paths_graph
from matplotlib import pyplot as plt
from indra.explanation import cycle_free_paths

num_nodes = 6
source = 0
target = num_nodes - 1

def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')

# Create a random graph
edge_prob = 0.5
randg = nx.random_graphs.erdos_renyi_graph(num_nodes, edge_prob, directed=True)
draw(randg, 'randg_%d_%s.pdf' % (num_nodes, edge_prob))

# ---------------------
# Use networkx to count all paths up to a certain length
nx_start = time.time()
nx_paths = []
try:
    nx_paths = list(nx.all_simple_paths(randg, source, target))
except nx.NetworkXNoPath:
    num_nx_paths = 0
nx_paths = set([tuple(p) for p in nx_paths])
nx_path_lengths = [len(p) for p in nx_paths]
nx_end = time.time()
nx_elapsed = nx_end - nx_start
print('NX time: %s' % nx_elapsed)
plt.close('all')
plt.ion()
bins = range(num_nodes+1)
nxh = plt.hist(nx_path_lengths, bins=bins, color='r', align='left')
# ---------------------

print("Starting PG")
pg_start = time.time()


(f_level, b_level) = paths_graph.get_reachable_sets(randg, source, target,
                            max_depth=num_nodes+1, signed=False)

sp_paths = []
for length in range(1, num_nodes):
    pg_raw = paths_graph.paths_graph(randg, source, target, length, f_level,
                                     b_level, signed=False)
    src = (0, source)
    tgt = (length, target)
    pg_0 = cycle_free_paths.PG_0(pg_raw, src, tgt)
    dic_PG = cycle_free_paths.PG(pg_0, src, tgt, length)
    G_cf, T = dic_PG[length - 1]
    P = cycle_free_paths.cf_sample_many_paths(src, tgt, G_cf, T, 10000)
    sp_paths += set(P)


#sp_lengths = [len(p) for p in sp_paths]
unique_paths = set(sp_paths)
up_lengths = [len(p) for p in unique_paths]

pg_end = time.time()
pg_elapsed = pg_end - pg_start
print('PG time: %s' % pg_elapsed)
sph = plt.hist(up_lengths, bins=bins, color='b', alpha=0.5,
               align='left')
#plt.figure()
#plt.hist(sp_lengths, bins=bins, color='g', align='left')
print(nx_paths)
print(unique_paths)
print(nx_paths == unique_paths)

"""
# Possible path lengths range from 1 (nodes neighbor each other) to num_nodes-1
# (nodes are arranged as a linear chain).
num_paths = 0
max_path_length = 10
for path_length in range(1, num_nodes):
    pg = paths_graph.paths_graph(randg, source, target, path_length, f_level,
                                 b_level, signed=False)
    num_cur_paths = 0
    if pg:
        paths = list(nx.all_simple_paths(
                                pg, (path_length, source), (0, target)))
        # Ignore cycles
        for path in paths:
            if len(set([n[1] for n in path])) == path_length + 1:
                num_cur_paths += 1
    num_paths += num_cur_paths
    print('%d paths of length %d' % (num_cur_paths, path_length))
print (pg_elapsed / nx_elapsed)

assert num_paths == num_nx_paths
"""

