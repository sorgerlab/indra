from indra.preassembler import grounding_mapper as gm
import sys
import logging
from indra.statements import Complex
import pickle
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.databases import biogrid_client as bg
import csv

logger = logging.getLogger('analyze_complexes')

# Filter out complexes from the statement list
def load_file(stmts_file):
    logger.info("Loading results...")
    with open(stmts_file) as f:
        results = pickle.load(f)
    return results


def get_biogrid_interactors(gene_list):
    pass

if __name__ == '__main__':

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

    # Get complexes
    complexes = [s for s in pa.unique_stmts if isinstance(s, Complex)]
    protein_complexes = [s for s in complexes
                           if all([True if 'HGNC' in ag.db_refs.keys()
                                        else False
                                        for ag in s.agent_list()])]

    genes = list(set([ag.db_refs['HGNC'] for stmt in protein_complexes
                                         for ag in stmt.members]))

    # Get complexes from BioGrid and combine duplicates
    bg_complexes = bg.get_statements(genes)
    pa_bg = Preassembler(hierarchies, bg_complexes)
    pa_bg.combine_duplicates()
    bg_complexes = pa_bg.unique_stmts

    reach_bg_matches = []
    for reach_stmt in protein_complexes:
        bg_matches = []
        for bg_stmt in bg_complexes:
            if bg_stmt.refinement_of(reach_stmt, hierarchies):
                bg_matches.append(bg_stmt.evidence[0].source_id)
                print 'BG %s refines/matches REACH %s' % (bg_stmt, reach_stmt)
        if bg_matches:
            reach_bg_matches.append((reach_stmt, bg_matches))

    not_matched = set(protein_complexes).difference(
                        set([t[0] for t in reach_bg_matches]))
    not_matched = sorted(not_matched, key=lambda x: len(x.evidence),
                         reverse=True)

    rows = []
    for stmt in not_matched:
        rows.append([stmt.members[0].name, stmt.members[1].name,
                     len(stmt.evidence)])
    with open('unmatched_complexes.tsv', 'w') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        csvwriter.writerows(rows)

