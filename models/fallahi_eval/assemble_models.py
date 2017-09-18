import os
from indra.util import _require_python3
import indra.tools.assemble_corpus as ac
from indra.tools.gene_network import GeneNetwork
from indra.tools.reading.submit_reading_pipeline_aws import \
    submit_run_reach, submit_combine, wait_for_complete
from assemble_pysb import *
import process_data

# This is the base name used for all files created/saved
basen = 'fallahi_eval'
# This is the base folder to read/write (potentially large) files from/to
# MODIFY ACCORDING TO YOUR OWN SETUP
based = os.path.join(os.environ['HOME'], 'data', basen)


def run_reading(pmid_fname):
    """Run reading on Amazon to produce a pickle of INDRA Statements

    This needs to be done only once, after that the resulting pickle file
    can be loaded.
    """
    # Submit reading
    job_list = submit_run_reach(basen, pmid_fname)
    # Wait for reading to complete
    reading_res = wait_for_complete(job_list)
    # Submid pickle combine job
    combine_res = submit_combine(basen, job_list)
    # Download the file and save


def build_prior(gene_names):
    """Build a corpus of prior Statements from PC and BEL."""
    gn = GeneNetwork(gene_names, basen)
    stmts = gn.get_statements(filter=False)
    ac.dump_statements(stmts, os.path.join(based, basen + '_prior.pkl'))


if __name__ == '__main__':
    data = process_data.read_rppa_data()
    gene_names = process_data.get_gene_names(data)

    # If generic assembly needs to be done (instead of just loading the result)
    # set this to True
    reassemble = False

    raw_stmts_file = os.path.join(based, basen + '_reach.pkl')
    pre_stmts_file = os.path.join(based, basen + '_reach_preassembled.pkl')

    if reassemble:
        # Load the raw statements
        stmts = ac.load_statements(raw_stmts_file)
        # Fix grounding and filter to grounded entities and for proteins,
        # filter to the human ones
        stmts = ac.map_grounding(stmts)
        stmts = ac.filter_grounded_only(stmts)
        stmts = ac.filter_human_only(stmts)
        # Combinatorially expand protein families
        stmts = ac.expand_families(stmts)
        # Apply a strict filter to statements based on the gene names
        stmts = ac.filter_gene_list(stmts, gene_names, 'all')
        # Fix errors in references to protein sequences
        stmts = ac.map_sequence(stmts)
        # Run preassembly and save result
        stmts = ac.run_preassembly(stmts, return_toplevel=False)
        ac.dump_statements(stmts, pre_stmts_file)

    # Load the preassembled statements
    stmts = ac.load_statements(pre_stmts_file)
    # Run assembly into a PySB model
    model = assemble_pysb(stmts, gene_names)
