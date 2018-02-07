from __future__ import absolute_import, print_function, unicode_literals, \
                       division
from builtins import dict, str
import sys
import subprocess
from argparse import ArgumentParser

if __name__ == '__main__':
    parser = ArgumentParser(
        description = 'This is a script designed to be used only with aws.'
        )
    parser.add_argument('pmid_list')
    parser.add_argument('tmp_dir')
    parser.add_argument('num_nodes', type=int)
    parser.add_argument('num_cores', type=int)
    parser.add_argument(
        '-r', '--readers',
        dest='readers',
        default='all',
        nargs='+',
        )
    args = parser.parse_args()

    READING_JOB_NAME = 'runreading'
    PROCESS_JOB_NAME = 'processreading'
    ASSEMBLE_JOB_NAME = 'assemblereading'

    # First, figure out how many papers there are by counting the number of
    # lines in the file.
    with open(args.pmid_list) as f:
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
    if num_pmids // args.num_nodes == num_pmids / args.num_nodes:
        node_chunk_size = int(num_pmids / args.num_nodes)
    else:
        node_chunk_size = int((num_pmids // args.num_nodes) + 1)
    # Limit jobs to a maximum chunk size to ensure regular checkpointing
    """
    MAX_NODE_CHUNK_SIZE = 2000
    if node_chunk_size > MAX_NODE_CHUNK_SIZE:
        node_chunk_size = MAX_NODE_CHUNK_SIZE
    """
    node_start_pts = range(0, num_pmids, node_chunk_size)

    for node_start_ix in node_start_pts:
        node_end_ix = node_start_ix + node_chunk_size

        cmd_list = [
            'qsub',
            '-b', 'y', '-V', '-cwd', '-N', READING_JOB_NAME, '-pe', 'orte',
            str(args.num_cores),
            'python',
            '-m', 'indra.tools.reading.pmid_reading.read_pmids',
            '-n', str(args.num_cores),
            '-s', str(node_start_ix),
            '-e', str(node_end_ix),
            args.tmp_dir,
            args.pmid_list,
            '-r'] + args.readers
        print(' '.join(cmd_list))
        subprocess.call(cmd_list)

    # Next, submit jobs for processing the REACH statements with INDRA. These
    # will be parallelized across cores rather than nodes, so we divide the
    # jobs accordingly. We make the start of these jobs contingent on the
    # completion of all of the REACH jobs.
    if num_pmids // (args.num_cores * args.num_nodes) == \
       num_pmids / (args.num_cores * args.num_nodes):
        core_chunk_size = int(num_pmids / (args.num_cores * args.num_nodes))
    else:
        core_chunk_size = int((num_pmids // (args.num_cores * args.num_nodes)) + 1)
    core_start_pts = range(0, num_pmids, core_chunk_size)

    for core_start_ix in core_start_pts:
        core_end_ix = core_start_ix + core_chunk_size
        cmd_list = [
            'qsub', '-b', 'y', '-V', '-cwd', '-hold_jid',
            READING_JOB_NAME, '-N', PROCESS_JOB_NAME, 'python', '-m',
            'indra.tools.reading.starcluster_reading.process_reach_from_s3', args.pmid_list,
            str(core_start_ix), str(core_end_ix)
            ]
        print(' '.join(cmd_list))
        subprocess.call(cmd_list)

    # Finally, we submit a job, contingent on the processing job, that
    # assembles all of the pickle files produced into a single file.
    cmd_list = ['qsub', '-b', 'y', '-V', '-cwd', '-hold_jid', PROCESS_JOB_NAME,
                '-N', ASSEMBLE_JOB_NAME, '-m', 'aes', '-M',
                'bachmanjohn@gmail.com', 'python', '-m',
                'indra.tools.reading.assemble_reading_stmts',
                'reach_stmts_*.pkl']
    print(' '.join(cmd_list))
    subprocess.call(cmd_list)
