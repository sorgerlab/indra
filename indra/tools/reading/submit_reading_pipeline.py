from __future__ import absolute_import, print_function, unicode_literals, \
                       division
from builtins import dict, str
import sys
import subprocess

if __name__ == '__main__':
    usage = 'Usage: %s pmid_list tmp_dir num_nodes num_cores_per_node' % \
             sys.argv[0]
    if len(sys.argv) != 5:
        print(usage)
        sys.exit()

    # The file containing the PMIDs to read
    pmid_list = sys.argv[1]
    # Path to temporary directory
    tmp_dir = sys.argv[2]
    # The number of nodes
    num_nodes = int(sys.argv[3])
    # Number of cores per node
    num_cores = int(sys.argv[4])

    REACH_JOB_NAME = 'runreach'
    PROCESS_JOB_NAME = 'processreach'
    ASSEMBLE_JOB_NAME = 'assemblereach'

    # First, figure out how many papers there are by counting the number of
    # lines in the file.
    with open(pmid_list) as f:
        for i, line in enumerate(f):
            pass
    num_pmids = i + 1
    if num_pmids == 0:
        print("No papers in PMID list.")
        sys.exit(1)

    # Next, submit jobs for offline reading with REACH. This will involve
    # running num_nodes jobs, each with num_cores assigned to each. The papers
    # are divided up into chunks such that all the papers will be processed in
    # a single round.
    if num_pmids // num_nodes == num_pmids / num_nodes:
        node_chunk_size = int(num_pmids / num_nodes)
    else:
        node_chunk_size = int((num_pmids // num_nodes) + 1)
    node_start_pts = range(0, num_pmids, node_chunk_size)

    for node_start_ix in node_start_pts:
        node_end_ix = node_start_ix + node_chunk_size

        cmd_list = ['qsub', '-b', 'y', '-V', '-cwd', '-N', REACH_JOB_NAME,
                    '-pe', 'orte', str(num_cores), 'python', '-m',
                    'indra.tools.reading.run_reach_on_pmids', pmid_list,
                    tmp_dir, str(num_cores), str(node_start_ix),
                    str(node_end_ix)]
        print(' '.join(cmd_list))
        subprocess.call(cmd_list)

    # Next, submit jobs for processing the REACH statements with INDRA. These
    # will be parallelized across cores rather than nodes, so we divide the
    # jobs accordingly. We make the start of these jobs contingent on the
    # completion of all of the REACH jobs.
    if num_pmids // (num_cores * num_nodes) == \
       num_pmids / (num_cores * num_nodes):
        core_chunk_size = int(num_pmids / (num_cores * num_nodes))
    else:
        core_chunk_size = int((num_pmids // (num_cores * num_nodes)) + 1)
    core_start_pts = range(0, num_pmids, core_chunk_size)

    for core_start_ix in core_start_pts:
        core_end_ix = core_start_ix + core_chunk_size
        cmd_list = ['qsub', '-b', 'y', '-V', '-cwd', '-hold_jid',
                    REACH_JOB_NAME, '-N', PROCESS_JOB_NAME, 'python', '-m',
                    'indra.tools.reading.process_reach_from_s3', pmid_list,
                    str(core_start_ix), str(core_end_ix)]
        print(' '.join(cmd_list))
        subprocess.call(cmd_list)

    # Finally, we submit a job, contingent on the processing job, that assembles
    # all of the pickle files produced into a single file.
    cmd_list = ['qsub', '-b', 'y', '-V', '-cwd', '-hold_jid', PROCESS_JOB_NAME,
                '-N', ASSEMBLE_JOB_NAME, '-m', 'aes', '-M',
                'bachmanjohn@gmail.com', 'python', '-m',
                'indra.tools.reading.assemble_reach_stmts',
                'reach_stmts_*.pkl']
    print(' '.join(cmd_list))
    subprocess.call(cmd_list)

