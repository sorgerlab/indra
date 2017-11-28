import sys
import random
import pickle
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
import pybel
from indra.util import _require_python3
from indra.util import plot_formatting as pf
from indra.tools import assemble_corpus as ac
from indra.explanation import paths_graph as pg
from indra.assemblers import PybelAssembler
from indra.util import read_unicode_csv, write_unicode_csv
from indra.assemblers.pysb_assembler import PysbPreassembler

def stmts_to_digraph(stmts):
    digraph = nx.MultiDiGraph()
    for stmt in stmts:
        agent_names = [a.name for a in stmt.agent_list() if a is not None]
        if len(agent_names) != 2:
            continue
        digraph.add_edge(agent_names[0], agent_names[1],
                         attr_dict={'uuid': stmt.uuid})
    return digraph


def get_pybel_node(pybel_model, gene_name):
    nodes = [node for node in pybel_model.nodes()
             if len(node) == 3 and \
                node[0] == pybel.constants.PROTEIN and \
                node[1] == 'HGNC' and \
                node[2] == gene_name]
    if len(nodes) == 0:
        return None
    elif len(nodes) == 1:
        return nodes[0]
    assert False


def stmts_to_pybel_graph(stmts):
    pba = PybelAssembler(stmts, name='INDRA/REACH Fallahi Eval Model',
                         description='Automatically assembled model.',
                         version='0.0.1')
    pba.make_model()
    return pba.to_signed_graph()


def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')


def filter_stmts(g, source, target, max_depth=6, signed=False,
                 target_polarity=None):
    f_level, b_level = pg.get_reachable_sets(g, source, target,
                                             max_depth=max_depth,
                                             signed=signed)
    stmt_uuids = set()
    stmt_nodes = set()
    stmt_uuid_nums = []
    stmt_node_nums = []
    # Iterate over various path lengths
    for length in range(1, max_depth+1):
        #print("Generating paths_graph for length %d" % length)
        this_pg = pg.paths_graph(g, source, target, length, f_level, b_level,
                                 signed=signed, target_polarity=target_polarity)
        if this_pg and len(this_pg) < 100:
            draw(this_pg, 'pybel_pg_%d.pdf' % length)

        # Get nodes for this length PG
        nodes_this_length = set([n[1] for n in this_pg])
        # Get stmt UUIDs for this length PG
        stmt_uuids_this_length = set()
        """
        for u, v in this_pg.edges():
            u_name, v_name = (u[1], v[1])
            multiedge_data = g.get_edge_data(u_name, v_name)
            for edge_data in multiedge_data.values():
                stmt_uuids_this_length.add((u_name, edge_data['uuid'], v_name))
        stmt_uuids |= stmt_uuids_this_length
        """
        # Get counts for this depth
        # Terminate the loop if we've saturated the number of edges
        #if stmt_nodes and stmt_nodes == nodes_this_length:
        #    break
        stmt_nodes |= nodes_this_length

        stmt_uuid_nums.append(len(stmt_uuids))
        stmt_node_nums.append(len(stmt_nodes))

    return (stmt_uuids, stmt_nodes, stmt_node_nums, stmt_uuid_nums)


def filter_stmts_all_paths(g, source, target):
    all_paths = nx.all_simple_paths(g, source, target)
    nodes = set()
    for path in all_paths:
        nodes |= set([n for n in path])
        print(len(nodes))


def filter_stmts_sampling(g, source, target, max_depth):
    f_level, b_level = pg.get_reachable_sets(g, source, target,
                                             max_depth=max_depth, signed=False)
    for length in range(1, max_depth+1):
        #print("Generating paths_graph for length %d" % length)
        this_pg = pg.paths_graph(g, source, target, length, f_level, b_level,
                                 signed=False)


def plot_results(g, stmt_uuids, stmt_nodes, stmt_node_nums, stmt_uuid_nums):
    norm_node_counts = np.array(stmt_node_nums) / len(g)
    norm_uuid_counts = np.array(stmt_uuid_nums) / len(g.edges())

    pf.set_fig_params()
    plt.ion()
    lengths = range(1, len(norm_uuid_counts)+1)
    plt.figure(figsize=(2, 2), dpi=150)
    plt.plot(lengths, norm_uuid_counts, color='orange', alpha=0.8,
             label='Statements')
    plt.plot(lengths, norm_node_counts, color='blue', alpha=0.8, label='Nodes')
    plt.legend(loc='upper left', fontsize=pf.fontsize, frameon=False)
    ax = plt.gca()
    pf.format_axis(ax)


if __name__ == '__main__':
    source = sys.argv[2]
    target = sys.argv[3]
    if len(sys.argv) > 4:
        max_depth = int(sys.argv[4])
    stmts = ac.load_statements(sys.argv[1])
    print(len(stmts))
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    stmts = [s for s in stmts if s.agent_list()[0]]
    print(len(stmts))
    from util import pkldump
    import ipdb; ipdb.set_trace()

    #ppa = PysbPreassembler(stmts)
    #ppa.replace_activities()
    #stmts = ppa.statements

    #g = stmts_to_digraph(stmts)
    g = stmts_to_pybel_graph(stmts)
    scc_lens = [len(s) for s in nx.strongly_connected_components(g)]
    scc_lens.sort(reverse=True)
    print("Largest strongly connected components: %s" % str(scc_lens[0:3]))
    # -- For PyBEL model --
    source = get_pybel_node(g, source)
    target = get_pybel_node(g, target)
    signed = True
    target_polarity = 0
    assert source and target

    for i in range(1):
        #source = random.choice(g.nodes())
        #target = random.choice(g.nodes())
        results = filter_stmts(g, source, target, max_depth=max_depth,
                               signed=signed, target_polarity=target_polarity)
        #filter_stmts_all_paths(g, source, target)
        plot_results(g, *results)
        #print(source)
        #print(target)
        #print("Final number of nodes: %d" % len(results[1]))
        #print("Final number of edges: %d" % len(results[0]))
