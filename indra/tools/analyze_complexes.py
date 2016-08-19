from indra.preassembler import grounding_mapper as gm
import sys
import logging
from indra.statements import Complex
import pickle
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies

logger = logging.getLogger('analyze_complexes')

# Filter out complexes from the statement list
def load_file(stmts_file):
    logger.info("Loading results...")
    with open(stmts_file) as f:
        results = pickle.load(f)
    return results


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

