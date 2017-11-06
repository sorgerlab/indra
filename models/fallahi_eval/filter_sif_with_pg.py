from indra.util import read_unicode_csv, write_unicode_csv
import networkx as nx
from indra.explanation import paths_graph as pg

sif_filename = 'output/fallahi_eval_preassembled_uuid_filtered.sif'

g = nx.MultiDiGraph()
for row in read_unicode_csv(sif_filename, delimiter=','):
    if len(row) != 3:
        continue
    source, uuid, target = row
    g.add_edge(source, target, attr_dict={'uuid': uuid})

source = 'BRAF'
target = 'JUN'

max_depth = 10
f_level, b_level = pg.get_reachable_sets(g, source, target, max_depth=max_depth,
                                         signed=False)
stmt_uuids = set()

for length in range(1, max_depth):
    print("Generating paths_graph for length %d" % length)
    this_pg = pg.paths_graph(g, source, target, length, f_level, b_level,
                             signed=False)

    stmt_uuids_this_length = set()
    for u, v in this_pg.edges():
        u_name, v_name = (u[1], v[1])
        multiedge_data = g.get_edge_data(u_name, v_name)
        for edge_data in multiedge_data.values():
            stmt_uuids_this_length.add((u_name, edge_data['uuid'], v_name))
    stmt_uuids |= stmt_uuids_this_length
    print("Paths of length %d: %d uuids" % (length, len(stmt_uuids_this_length)))

write_unicode_csv('fallahi_eval_paths_graph_BRAF_JUN_max%d' % max_depth,
                  list(stmt_uuids))
