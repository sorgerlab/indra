from __future__ import absolute_import, print_function, unicode_literals, \
                       division
from builtins import dict, str
import sys
import subprocess
import numpy as np

if __name__ == '__main__':
    usage = 'Usage: %s pmid_list num_jobs'
    if len(sys.argv) != 3:
        print(usage)
        sys.exit()

    pmid_list = sys.argv[1]
    num_jobs = int(sys.argv[2])

    with open(pmid_list) as f:
        pmids = [line.strip('\n') for line in f.readlines()]

    ids_per_job = len(pmids) / float(num_jobs)
    start_indices = np.arange(0, len(pmids), ids_per_job)
    for ix, start_ix in enumerate(start_indices):
        start_ix = int(start_ix)
        if ix + 1 < num_jobs:
            end_ix = int(start_indices[ix + 1])
        else:
            end_ix = len(pmids)
        out_filename = '%s_%d_%d.out' % (pmid_list, start_ix, end_ix)
        cmd_list = ['bsub', '-q', 'short', '-W', '12:00', '-N', '-o',
                    out_filename, 'python', '-m',
                    'indra.tools.reading.upload_content_to_s3',
                    pmid_list, str(start_ix), str(end_ix), 'force_fulltext']
        print(' '.join(cmd_list))
        subprocess.call(cmd_list)

