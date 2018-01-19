from indra.util import _require_python3, write_unicode_csv
import indra.tools.assemble_corpus as ac
#from indra.tools.reading.submit_reading_pipeline import \
#    submit_run_reach, submit_combine, wait_for_complete
import assemble_pysb
import process_data
from util import *


def run_reading(pmid_fname):
    """Run reading on Amazon to produce a pickle of INDRA Statements

    This needs to be done only once, after that the resulting pickle file
    can be loaded.
    """
    # Submit reading
    job_list = submit_run_reach(basen, pmid_fname)
    # Wait for reading to complete
    reading_res = wait_for_complete(job_list)
    # Submit pickle combine job
    combine_res = submit_combine(basen, job_list)
    # Download the file and save


def get_stmt_sif(stmts, fname):
    rows = []
    for stmt in stmts:
        agent_names = [a.name for a in stmt.agent_list() if a is not None]
        if len(agent_names) != 2:
            continue
        rows.append((agent_names[0], stmt.uuid, agent_names[1]))
    write_unicode_csv(fname, rows)


if __name__ == '__main__':
    # Load the data and get the gene names
    data = process_data.read_rppa_data()
    gene_names = process_data.get_gene_names(data)

    # If generic assembly needs to be done (instead of just loading the result)
    # set this to True
    reassemble = False

    # The file in which the preassembled statements will be saved
    pre_stmts_file = prefixed_pkl('preassembled')
    if reassemble:
        # Load various files that were previously produced
        sources = ['indradb', 'trips', 'bel', 'biopax', 'phosphosite', 'r3',
                   'sparser']
        stmts = []
        for source in sources:
            stmts += ac.load_statements(prefixed_pkl(source))
        stmts = ac.filter_no_hypothesis(stmts)
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
    assemble_pysb.assemble_pysb(stmts, gene_names, contextualize=True)
