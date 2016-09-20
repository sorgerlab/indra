#import rdflib
#from rdflib import Namespace, Literal
import sys
import logging
import pickle
from indra.preassembler import grounding_mapper as gm
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
import urllib2, urllib
import cStringIO
import sys
import csv
from indra.tools import plot_formatting as pf
from matplotlib import pyplot as plt
import numpy as np

logger = logging.getLogger('analyze_biological_processes')


def load_file(stmts_file):
    logger.info("Loading results...")
    with open(stmts_file) as f:
        results = pickle.load(f)
    return results


# Get statements involving biological processes
def go_protein_pair(stmt):
    go = None
    protein = None
    for ag in stmt.agent_list():
        if ag is None:
            continue
        grounding = ag.db_refs.keys()
        if 'HGNC' in grounding:
            protein = ag.db_refs.get('HGNC')
        elif 'GO' in grounding:
            go = ag.db_refs.get('GO')
            bp_name = ag.name
    if go is not None and protein is not None:
        return (bp_name, go, protein)
    else:
        return (None, None, None)


def get_genes_for_go_id(goid):
    quickgo_url = 'https://www.ebi.ac.uk/QuickGO/GAnnotation'
    params = {'goid': goid, 'format':'tsv', 'db':'UniProtKB', 'tax':'9606',
              'col':'proteinSymbol'
           }
    try:
        res = urllib2.urlopen(quickgo_url, data=urllib.urlencode(params))
    except urllib2.HTTPError:
        logging.error('Could not retrieve proteins associated with GO ID %s'
                      % goid)
        return None
    tsv_str = cStringIO.StringIO(res.read())
    tsv_reader = csv.reader(tsv_str, delimiter='\t')
    genes = set([])
    for row in tsv_reader:
        genes.add(row[0])
    return list(genes)


def plot_stmt_counts(go_stmt_map, plot_filename, figsize=(3, 3)):
    # Put together counts for a figure
    pf.set_fig_params()
    fig = plt.figure(figsize=(8, 4), dpi=200)
    ax = fig.gca()
    counts = []
    for go_id in go_stmt_map.keys():
        counts.append((go_stmt_map[go_id]['names'][0],
                       len(go_stmt_map[go_id]['in_go']),
                       len(go_stmt_map[go_id]['not_in_go'])))
    counts.sort(key=lambda x: x[1] + x[2], reverse=True)
    indices = np.arange(len(counts))
    width = 0.8
    in_go = [c[1] for c in counts]
    labels = [c[0] for c in counts]
    not_in_go = [c[2] for c in counts]
    p1 = ax.bar(indices, in_go, width, color='b')
    p2 = ax.bar(indices, not_in_go, width, color='r', bottom=in_go)
    ax.set_ylabel('No. of stmts')
    plt.xlim([-0.2, len(counts)])
    plt.xticks(indices + width/2., labels, rotation='vertical')
    plt.legend((p1[0], p2[0]), ('In GO', 'Not in GO'), frameon=False,
               fontsize=7)
    plt.subplots_adjust(left=0.08, bottom=0.44, top=0.85, right=0.97)
    pf.format_axis(ax)
    fig.savefig(plot_filename)

if __name__ == '__main__':

    #with open('go_stmt_map.pkl') as f:
    #    import pickle
    #    go_stmt_map = pickle.load(f)
    #plot_stmt_counts(go_stmt_map, 'go_stmts.pdf')
    #sys.exit()

    # Load the statements
    if len(sys.argv) < 2:
        print "Usage: %s reach_stmts_file" % sys.argv[0]
        sys.exit()
    results = load_file(sys.argv[1])

    all_stmts = [stmt for paper_stmts in results.values()
                      for stmt in paper_stmts]

    # Map grounding
    logger.info('Mapping grounding...')
    gmap = gm.GroundingMapper(gm.default_grounding_map)
    map_stmts = gmap.map_agents(all_stmts)
    map_stmts = gmap.rename_agents(map_stmts)

    # Combine duplicates
    logger.info('Removing duplicates...')
    pa = Preassembler(hierarchies, map_stmts)
    pa.combine_duplicates()

    # Map GO IDs to genes and associated statements
    logger.info('Building map from GO IDs to stmts')
    go_protein_map = {}
    go_name_map = {}
    for stmt in pa.unique_stmts:
        (bp_name, go, hgnc) = go_protein_pair(stmt)
        if bp_name is None and go is None and hgnc is None:
            continue
        go_prot_list = go_protein_map.get(go, [])
        go_prot_list.append((hgnc, stmt))
        go_protein_map[go] = go_prot_list
        go_name_set = go_name_map.get(go, set([]))
        go_name_set.add(bp_name)
        go_name_map[go] = go_name_set

    # Iterate over all of the GO IDs and compare the annotated genes in GO
    # to the ones from the given statements
    go_stmt_map = {}
    for ix, go_id in enumerate(go_protein_map.keys()):
        logger.info('Getting genes for %s (%s) from GO (%d of %d)' %
                    (go_id, ','.join(list(go_name_map[go_id])),
                     ix+1, len(go_protein_map.keys())))
        genes_from_go = get_genes_for_go_id(go_id)
        prot_stmt_list = go_protein_map[go_id]
        in_go = []
        not_in_go = []
        for (prot, stmt) in prot_stmt_list:
            if prot in genes_from_go:
                in_go.append(stmt)
            else:
                not_in_go.append(stmt)
        go_stmt_map[go_id] = {'names': list(go_name_map[go_id]),
                              'in_go': in_go, 'not_in_go': not_in_go}

    with open('go_stmt_map.pkl', 'w') as f:
        pickle.dump(go_stmt_map, f)

    plot_stmt_counts(go_stmt_map, 'go_stmts.pdf')
