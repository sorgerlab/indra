from indra.util import read_unicode_csv, write_unicode_csv
from indra.util import plot_formatting as pf
import networkx as nx
from indra.explanation import paths_graph as pg
from matplotlib import pyplot as plt
from indra.util import _require_python3
import numpy as np
import random
import pickle

sif_filename = 'output/fallahi_eval_preassembled_uuid_filtered.sif'

# Load directed and undirected versions of the SIF
ug = nx.Graph()
g = nx.MultiDiGraph()
for row in read_unicode_csv(sif_filename, delimiter=','):
    if len(row) != 3:
        continue
    source, uuid, target = row
    g.add_edge(source, target, attr_dict={'uuid': uuid})
    ug.add_edge(source, target)

scc_sizes = [len(scc) for scc in nx.strongly_connected_components(g)]

for i in range(1):
    #source = random.choice(ug.nodes())
    #target = random.choice(ug.nodes())
    source = 'BRAF'
    target = 'JUN'

    max_depth = 8
    total_uuids = len(g.edges())
    total_nodes = len(g)

    f_level, b_level = pg.get_reachable_sets(g, source, target,
                                             max_depth=max_depth, signed=False)

    stmt_uuids = set()
    stmt_nodes = set()
    stmt_uuid_nums = []
    stmt_node_nums = []

    for length in range(1, max_depth+1):
        print("Generating paths_graph for length %d" % length)
        this_pg = pg.paths_graph(g, source, target, length, f_level, b_level,
                                 signed=False)

        # Get nodes for this length PG
        nodes_this_length = set([n[1] for n in this_pg])
        # Get stmt UUIDs for this length PG
        stmt_uuids_this_length = set()
        for u, v in this_pg.edges():
            u_name, v_name = (u[1], v[1])
            multiedge_data = g.get_edge_data(u_name, v_name)
            for edge_data in multiedge_data.values():
                stmt_uuids_this_length.add((u_name, edge_data['uuid'], v_name))
        stmt_uuids |= stmt_uuids_this_length
        stmt_nodes |= nodes_this_length
        # Get counts for this depth
        if stmt_uuid_nums and stmt_uuid_nums[-1] != 0 and \
                len(stmt_uuids) == stmt_uuid_nums[-1]:
            break
        stmt_uuid_nums.append(len(stmt_uuids))
        stmt_node_nums.append(len(stmt_nodes))
        print("Paths of length %d: %d uuids" %
              (length, len(stmt_uuids_this_length)))

    write_unicode_csv('fallahi_eval_paths_graph_BRAF_JUN_max%d' % max_depth,
                      list(stmt_uuids))

    norm_node_counts = np.array(stmt_node_nums) / total_nodes
    norm_uuid_counts = np.array(stmt_uuid_nums) / total_uuids

    pf.set_fig_params()

    plt.ion()
    lengths = range(len(norm_uuid_counts))
    plt.figure(figsize=(2, 2), dpi=150)
    plt.plot(lengths, norm_uuid_counts, color='orange', alpha=0.8,
             label='Statements')
    plt.plot(lengths, norm_node_counts, color='blue', alpha=0.8, label='Nodes')
    plt.legend(loc='upper left', fontsize=pf.fontsize, frameon=False)
    ax = plt.gca()
    pf.format_axis(ax)


