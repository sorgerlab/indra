import os
from indra.util import _require_python3
import indra.tools.assemble_corpus as ac
from indra.tools.reading.submit_reading_pipeline_aws import \
    submit_run_reach, wait_for_complete
from assemble_pysb import *
import process_data

def run_reading(pmid_fname):
    # Submit reading
    job_list = submit_run_reach('fallahi_eval', pmid_fname)
    # Wait for reading to complete
    reading_res = wait_for_complete(job_list)


if __name__ == '__main__':
    data = process_data.read_rppa_data()
    gene_names = process_data.get_gene_names(data)

    reassemble = False

    raw_stmts_file = '%s/data/fallahi_reach.pkl' % os.environ['HOME']
    pre_stmts_file = '%s/data/fallahi_reach_preassembled.pkl' % \
                        os.environ['HOME']

    if reassemble:
        stmts = ac.load_statements(raw_stmts_file)
        stmts = ac.map_grounding(stmts)
        stmts = ac.filter_grounded_only(stmts)
        stmts = ac.filter_human_only(stmts)
        stmts = ac.expand_families(stmts)
        stmts = ac.filter_gene_list(stmts, gene_names, 'all')
        stmts = ac.map_sequence(stmts)
        stmts = ac.run_preassembly(stmts, return_toplevel=False)
        ac.dump_statements(stmts, pre_stmts_file)

    stmts = ac.load_statements(pre_stmts_file)
    model = assemble_pysb(stmts, gene_names)
