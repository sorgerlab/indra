from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import pickle
import logging
from indra.statements import Complex
from indra.preassembler import Preassembler
from indra.preassembler import grounding_mapper as gm
from indra.preassembler.hierarchy_manager import hierarchies
from indra.databases import biogrid_client as bg
from indra.databases import hgnc_client, uniprot_client
from indra.util import write_unicode_csv

logger = logging.getLogger('complexes')

# Filter out complexes from the statement list
def load_file(stmts_file):
    logger.info("Loading results...")
    with open(stmts_file, 'rb') as f:
        results = pickle.load(f)
    return results

def analyze(filename):
    results = load_file(filename)

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
    # Get HGNC grounding
    protein_complexes = [s for s in complexes
                           if all([True if 'HGNC' in ag.db_refs.keys()
                                        else False
                                        for ag in s.agent_list()])]

    logger.info('Mapping gene IDs to gene symbols')
    gene_ids = list(set([ag.db_refs['HGNC'] for stmt in protein_complexes
                                            for ag in stmt.members]))
    genes = [hgnc_client.get_hgnc_name(id) for id in gene_ids]

    # Get complexes from BioGrid and combine duplicates
    num_genes_per_query = 50
    start_indices = range(0, len(genes), num_genes_per_query)
    end_indices = [i + num_genes_per_query
                   if i + num_genes_per_query < len(genes) else len(genes)
                   for i in start_indices]
    bg_complexes = []
    for i in range(len(start_indices)):
        logger.info("Querying biogrid for %s" %
                    str(genes[start_indices[i]:end_indices[i]]))
        bg_complexes += (bg.get_statements(
                                genes[start_indices[i]:end_indices[i]]))

    # Filter out Biogrid statements not involving genes in the gene list
    # (this will make duplicate removal more efficient
    bg_filt = []
    for stmt in bg_complexes:
        if stmt.members[0].name in genes and \
           stmt.members[1].name in genes:
            bg_filt.append(stmt)
    # Might as well free up some memory
    del bg_complexes

    logger.info("Combining duplicates with biogrid...")
    pa = Preassembler(hierarchies, bg_filt + protein_complexes)
    pa.combine_duplicates()

    indra_only = []
    bg_only = []
    indra_and_bg = []
    for stmt in pa.unique_stmts:
        evidence_source_list = set([])
        for e in stmt.evidence:
            evidence_source_list.add(e.source_api)
        if 'reach' in evidence_source_list and \
           'biogrid' in evidence_source_list:
            indra_and_bg.append(stmt)
        elif 'reach' in evidence_source_list and \
             'biogrid' not in evidence_source_list:
            indra_only.append(stmt)
        elif 'reach' not in evidence_source_list and \
             'biogrid' in evidence_source_list:
            bg_only.append(stmt)

    rows = []
    for stmt in indra_only:
        rows.append([stmt.members[0].name, stmt.members[1].name,
                     str(len(stmt.evidence))])
    write_unicode_csv('unmatched_complexes.tsv', rows, delimiter='\t')

    return {'indra_only': indra_only,
            'bg_only': bg_only,
            'indra_and_bg': indra_and_bg}

