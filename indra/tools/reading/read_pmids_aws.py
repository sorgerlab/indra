from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.tools.reading.read_pmids import READER_DICT

DOC = \
"""
This script is intended to be run on an Amazon ECS container, so information
for the job either needs to be provided in environment variables (e.g., the
REACH version and path) or loaded from S3 (e.g., the list of PMIDs).
"""

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description=DOC
        )
    parser.add_argument(
        dest='basename',
        help='The name of this run.'
        )
    parser.add_argument(
        dest='out_dir',
        help='The name of the temporary output directory'
        )
    parser.add_argument(
        dest='num_cores',
        help='Select the number of cores on which to run.',
        type=int
        )
    parser.add_argument(
        dest='start_index',
        help='Select the index of the first pmid in the list to read.',
        type=int
        )
    parser.add_argument(
        dest='end_index',
        help='Select the index of the last pmid in the list to read.',
        type=int
        )
    parser.add_argument(
        '-r', '--readers',
        dest='readers',
        default='all',
        choices=list(READER_DICT.keys()) + ['all'],
        nargs='+',
        help='Choose which reader(s) to use.'
        )
    args = parser.parse_args()
    from indra.tools.reading import read_pmids as read
    import boto3
    import botocore
    import os
    import pickle
    import logging
    import sys

    logger = logging.getLogger('read_pmids_aws')

    # Setting default force read/fulltext parameters
    force_read = False
    force_fulltext = False

    client = boto3.client('s3')
    bucket_name = 'bigmech'
    pmid_list_key = 'reading_results/%s/pmids' % args.basename
    path_to_reach = os.environ.get('REACH_JAR_PATH')
    reach_version = os.environ.get('REACH_VERSION')
    if path_to_reach is None or reach_version is None:
        print('REACH_JAR_PATH and/or REACH_VERSION not defined, exiting.')
        sys.exit(1)

    try:
        pmid_list_obj = client.get_object(
            Bucket=bucket_name,
            Key=pmid_list_key
            )
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == 'NoSuchKey':
            logger.info('Could not find PMID list file at %s, exiting' %
                        pmid_list_key)
            sys.exit(1)
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    # Get the content from the object
    pmid_list_str = pmid_list_obj['Body'].read().decode('utf8').strip()
    pmid_list = [line.strip() for line in pmid_list_str.split('\n')]

    # Handle the all option.
    if 'all' in args.readers:
        readers = list(READER_DICT.keys())
    else:
        readers = args.readers[:]

    # Run the reading pipelines
    stmts = {}
    content_types = {}
    for reader, run_reader in read.READER_DICT.items():
        if reader not in readers:
            continue

        some_stmts, some_content_types = run_reader(
                pmid_list,
                args.out_dir,
                args.num_cores,
                args.start_index,
                args.end_index,
                force_read,
                force_fulltext,
                cleanup=False,
                verbose=True
                )
        stmts[reader] = some_stmts
        content_types[reader] = some_content_types

    # Pickle the content types to S3
    for reader, some_stmts in stmts.items():
        N_papers = len(stmts[reader])
        ct_key_name = 'reading_results/%s/%s/content_types/%d_%d.pkl' % \
                      (args.basename, reader, args.start_index, args.end_index)
        logger.info("Saving content types for %d papers to %s" %
                    (N_papers, ct_key_name))
        ct_bytes = pickle.dumps(content_types[reader])
        client.put_object(Key=ct_key_name, Body=ct_bytes, Bucket=bucket_name)
        # Pickle the statements to a bytestring
        pickle_key_name = 'reading_results/%s/%s/stmts/%d_%d.pkl' % \
                          (args.basename, reader, args.start_index, args.end_index)
        logger.info("Saving stmts from %d papers to %s" %
                    (N_papers, pickle_key_name))
        stmts_bytes = pickle.dumps(stmts[reader])
        client.put_object(Key=pickle_key_name, Body=stmts_bytes,
                          Bucket=bucket_name)
