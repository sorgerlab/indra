from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import pickle
import logging
from indra.sources import reach
from indra.literature import s3_client

if __name__ == '__main__':
    usage = "Usage: %s pmid_list start_index end_index" % sys.argv[0]
    if len(sys.argv) < 4:
        print(usage)
        sys.exit()

    # Logger
    logger = logging.getLogger('processreach')

    pmid_list_file = sys.argv[1]
    start_ix = int(sys.argv[2])
    end_ix = int(sys.argv[3])

    # Load the list of PMIDs from the given file
    with open(pmid_list_file) as f:
        pmid_list = [line.strip('\n') for line in f.readlines()]
    if end_ix > len(pmid_list):
        end_ix = len(pmid_list)

    stmts = {}
    for ix, pmid in enumerate(pmid_list[start_ix:end_ix]):
        reach_json_str = s3_client.get_reader_json_str('reach', pmid)
        # Logging message will have been produced by get_reach_output
        if reach_json_str is None:
            continue
        # Run the REACH processor on the JSON
        try:
            logger.info('%d: Processing %s' % (ix, pmid))
            reach_proc = reach.process_json_str(reach_json_str, citation=pmid)
        # If there's a problem, skip it
        except Exception as e:
            print("Exception processing %s" % pmid)
            print(e)
            continue

        stmts[pmid] = reach_proc.statements

    with open('reach_stmts_%d_%d.pkl' % (start_ix, end_ix), 'wb') as f:
        pickle.dump(stmts, f, protocol=2)
