from indra.util import _require_python3
import os
import sys
from random import shuffle
from indra.sources import biopax
import indra.tools.assemble_corpus as ac
from indra.assemblers import CxAssembler
from indra.literature.pubmed_client import get_ids_for_gene


base_network = 'PAR1-mediated thrombin signaling events.owl'
base_path = os.path.join(os.environ['HOME'], 'data', 'nci_pid', base_network)

def process_owl(path):
    bp = biopax.process_owl(path)
    return bp.statements

def get_gene_names(stmts):
    names = set()
    for stmt in stmts:
        for agent in stmt.agent_list():
            if agent is not None and agent.db_refs.get('HGNC'):
                names.add(agent.name)
    return sorted(list(names))

def get_pmids(gene_names):
    pmids = []
    for gene_name in gene_names:
        pm = get_ids_for_gene(gene_name)
        pmids += pm
        print('%s: %d PMIDs' % (gene_name, len(pm)))
    return pmids

def save_pmids_for_reading(pmids, fname):
    shuffle(pmids)
    with open(fname, 'wt') as fh:
        for pmid in pmids:
            fh.write('%s\n' % pmid)

def get_reach_output(path):
    stmts = ac.load_statements(path)
    return stmts

def run_assembly(stmts, save_file):
    stmts = ac.map_grounding(stmts)
    stmts = ac.filter_grounded_only(stmts)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.expand_families(stmts)
    stmts = ac.filter_gene_list(stmts, gene_names, 'one')
    stmts = ac.map_sequence(stmts)
    stmts = ac.run_preassembly(stmts, return_toplevel=False)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_enzyme_kinase(stmts)
    ac.dump_statements(stmts, save_file)
    return stmts

def assemble_cx(stmts, save_file):
    cxa = CxAssembler(stmts)
    cxa.make_model(add_indra_json=False)
    cxa.save_model(save_file)
    return cxa

if __name__ == '__main__':
    biopax_stmts = process_owl(base_path)
    gene_names = get_gene_names(biopax_stmts)
    #pmids = get_pmids(gene_names)
    #save_pmids_for_reading(pmids, 'gpcr_pmids.txt')
    reach_stmts = get_reach_output('reach_output.pkl')
    stmts = biopax_stmts + reach_stmts
    stmts = run_assembly(stmts, 'final_stmts.pkl')
    cxa = assemble_cx(stmts, 'gpcr.cx')
