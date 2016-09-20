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
    if go is not None and protein is not None:
        return (go, protein)
    else:
        return (None, None)

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

def get_all_go_ids(stmt_list):
    go_ids = set([])
    for stmt in stmt_list:
        for ag in stmt.agent_list():
            if ag is not None and ag.db_refs.get('GO') is not None:
                go_ids.add(ag.db_refs.get('GO'))
    return list(go_ids)

if __name__ == '__main__':

    #gids = get_genes_for_go_id('GO:0006954')

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

    import ipdb;
    go_protein_map = {}
    for stmt in pa.unique_stmts:
        (go, hgnc) = go_protein_pair(stmt)
        if go is None and hgnc is None:
            continue
        go_prot_set = go_protein_map.get(go, set([]))
        go_prot_set.add(hgnc)
        go_protein_map[go] = go_prot_set

    # bp_stmts = [s for s in pa.unique_stmts if go_protein_stmt(s)]

    #genes = list(set([ag.db_refs.get('HGNC') for stmt in bp_stmts
    #                                         for ag in stmt.agent_list()
    #                                         if ag.db_refs.get('HGNC')]))
