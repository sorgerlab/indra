from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra import trips
from assembly_eval import have_file, run_assembly
import logging

if __name__ == '__main__':
    # Quiet the requests logging
    logging.getLogger('requests').setLevel(logging.ERROR)
    logging.getLogger('urllib3').setLevel(logging.ERROR)

    pmc_ids = [s.strip() for s in open('pmcids.txt', 'rt').readlines()]
    # Use the existing EKB extractions.
    for pmcid in pmc_ids:
        folder = 'trips'
        prefix = folder + '/' + pmcid
        print('Processing %s...' % pmcid)
        with open(prefix + '.ekb', 'r') as f:
            tp = trips.process_xml(f.read())
        # PMIDs from TRIPS need to be set here because it propagates
        # the PMCID by default
        run_assembly(tp.statements, folder, pmcid)
