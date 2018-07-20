"""This script is intended to be run on an Amazon ECS container, so information
for the job either needs to be provided in environment variables (e.g., the
REACH version and path) or loaded from S3 (e.g., the list of PMIDs).
"""
from __future__ import absolute_import, print_function, unicode_literals

import pickle
from builtins import dict, str
import boto3
import botocore
import logging
import sys
import os
import random
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from datetime import datetime
from indra.tools.reading.db_reading.read_db import produce_readings, \
    produce_statements, get_id_dict
from indra.tools.reading.readers import get_readers


def plot_hist(agged, agg_over, data, s3, s3_base, bucket_name):
    fig = plt.figure()
    plt.hist(data, bins=np.arange(len(data)))
    plt.xlabel('Number of %s for %s' % (agged, agg_over))
    plt.ylabel('Number of %s with a given number of %s' % (agg_over, agged))
    fname = '%s_per_%s.png' % (agged, agg_over)
    fig.savefig(fname)
    with open(fname, 'rb') as f:
        s3_key = s3_base + fname
        s3.put_object(Key=s3_key, Body=f.read(), Bucket=bucket_name)
    return


def report_statistics(reading_outputs, stmt_outputs, starts, ends,
                      s3_log_prefix, s3, bucket_name):
    s3_prefix = s3_log_prefix + 'statisics/'
    starts['stats'] = datetime.now()
    data_dict = {}
    hist_dict = {}
    data_dict['Total readings'] = len(reading_outputs)
    reading_stmts = [(rd.reading_id, rd.tcid, rd.reader, rd.get_statements())
                     for rd in reading_outputs]
    data_dict['Content processed'] = len({t[1] for t in reading_stmts})
    data_dict['Statements produced'] = len(stmt_outputs)

    readings_with_stmts = []
    readings_with_no_stmts = []
    for t in reading_stmts:
        if t[-1]:
            readings_with_stmts.append(t)
        else:
            readings_with_no_stmts.append(t)

    data_dict['Readings with 0 statements'] = len(readings_with_no_stmts)

    # Do a bunch of aggregation
    tc_rd_dict = {}
    tc_stmt_dict = {}
    rd_stmt_dict = {}
    reader_stmts = {}
    reader_tcids = {}
    reader_rids = {}
    for rid, tcid, reader, stmts in readings_with_stmts:
        # Handle things keyed by tcid
        if tcid not in tc_rd_dict.keys():
            tc_rd_dict[tcid] = {rid}
            tc_stmt_dict[tcid] = set(stmts)
        else:
            tc_rd_dict[tcid].add(rid)
            tc_stmt_dict[tcid] |= set(stmts)

        # Handle things keyed by rid
        if rid not in rd_stmt_dict.keys():
            rd_stmt_dict[rid] = set(stmts)
        else:
            rd_stmt_dict[rid] |= set(stmts)  # this shouldn't really happen.

        # Handle things keyed by reader
        if reader not in reader_stmts.keys():
            reader_stmts[reader] = set(stmts)
            reader_tcids[reader] = {tcid}
            reader_rids[reader] = {rid}
        else:
            reader_stmts[reader] |= set(stmts)
            reader_tcids[reader].add(tcid)
            reader_rids[reader].add(rid)

    for rid, tcid, reader, _ in readings_with_no_stmts:
        # Handle things keyed by tcid
        if tcid not in tc_rd_dict.keys():
            tc_rd_dict[tcid] = {rid}
        else:
            tc_rd_dict[tcid].add(rid)

        # Handle things keyed by reader
        if reader not in reader_stmts.keys():
            reader_tcids[reader] = {tcid}
            reader_rids[reader] = {rid}
        else:
            reader_tcids[reader].add(tcid)
            reader_rids[reader].add(rid)

    # Produce some numpy count arrays.
    hist_dict[('readings', 'text content')] = \
        np.array([len(rid_set) for rid_set in tc_rd_dict.values()])
    hist_dict[('statements', 'text_content')] = \
        np.array([len(stmts) for stmts in tc_stmt_dict.values()])
    hist_dict[('statements', 'readings')] = \
        np.array([len(stmts) for stmts in rd_stmt_dict.values()])
    hist_dict[('statements', 'readers')] = \
        np.array([len(stmts) for stmts in reader_stmts.values()])
    hist_dict[('text content', 'reader')] = \
        np.array([len(tcid_set) for tcid_set in reader_tcids.values()])
    hist_dict[('readings', 'reader')] = \
        np.array([len(rid_set) for rid_set in reader_rids.values()])

    # Produce the histograms
    for (agged, agg_over), data in hist_dict.items():
        plot_hist(agged, agg_over, data, s3, s3_prefix, bucket_name)
        label = '%s per %s' % (agged, agg_over)
        data_dict[label.capitalize()] = {'mean': data.mean(), 'std': data.std(),
                                         'median': np.median(data)}

    # Produce the summary report
    text_report_str = ''
    top_labels = ['Total readings', 'Content processed', 'Statements produced']
    for label in top_labels:
        text_report_str += '%s: %d\n' % (label, data_dict[label])

    for label, data in data_dict.items():
        if label in top_labels:
            continue
        if isinstance(data, dict):
            text_report_str += '%s:\n' % label
            text_report_str += '\n'.join(['\t%s: %d' % (k, v)
                                          for k, v in data.items()])
            text_report_str += '\n'
        else:
            text_report_str += '%s: %d\n' % (label, data)

    s3.put_object(Key=s3_prefix + 'summary.txt', Body=text_report_str,
                  Bucket=bucket_name)
    s3.put_object(Key=s3_prefix + 'hist_data.pkl', Body=pickle.dumps(hist_dict),
                  Bucket=bucket_name)
    s3.put_object(Key=s3_prefix + 'sum_data.pkl', Body=pickle.dumps(data_dict),
                  Bucket=bucket_name)
    ends['stats'] = datetime.now()

    # Report on the timing
    timing_str = ''
    for step in ['reading', 'statement production', 'stats']:
        time_taken = ends[step] - starts[step]
        timing_str += ('%22s: start: %s, end: %s, duration: %s\n'
                       % (step, starts[step], ends[step], time_taken))

    s3.put_object(Key=s3_prefix + 'time.txt', Body=timing_str, Bucket=bucket_name)
    return


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description=__doc__
        )
    parser.add_argument(
        dest='basename',
        help='The name of this run.'
        )
    parser.add_argument(
        dest='job_name',
        help='The name of this job.'
        )
    parser.add_argument(
        dest='out_dir',
        help='The name of the temporary output directory'
        )
    parser.add_argument(
        dest='read_mode',
        choices=['all', 'unread', 'none'],
        help=("Set the reading mode. If 'all', read everything, if "
              "'unread', only read content that does not have pre-existing "
              "readings of the same reader and version, if 'none', only "
              "use pre-existing readings. Default is 'unread'.")
        )
    parser.add_argument(
        dest='stmt_mode',
        choices=['all', 'unread', 'none'],
        help=("Choose which readings should produce statements. If 'all', all "
              "readings that are produced or retrieved will be used to produce "
              "statements. If 'unread', only produce statements from "
              "previously unread content. If 'none', do not produce any "
              "statements (only readings will be produced).")
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
        choices=['reach', 'sparser'],
        nargs='+',
        help='Choose which reader(s) to use.'
        )
    parser.add_argument(
        '--force_fulltext',
        action='store_true',
        help='Require that content be fulltext, skip anything that isn\'t.'
        )
    parser.add_argument(
        '--use_best_fulltext',
        action='store_true',
        help='Read only the "best" fulltext available for a given id.'
        )
    parser.add_argument(
        '--test',
        action='store_true',
        help="Use the test database."
        )
    args = parser.parse_args()

    logger = logging.getLogger('read_db_aws')
    logger.setLevel(logging.DEBUG)

    client = boto3.client('s3')
    bucket_name = 'bigmech'
    id_list_key = 'reading_inputs/%s/id_list' % args.basename
    readers = [reader_class(args.basename, args.num_cores)
               for reader_class in get_readers()
               if reader_class.name.lower() in args.readers]

    try:
        id_list_obj = client.get_object(
            Bucket=bucket_name,
            Key=id_list_key
            )
    # Handle a missing object gracefully
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == 'NoSuchKey':
            logger.info('Could not find PMID list file at %s, exiting' %
                        id_list_key)
            sys.exit(1)
        # If there was some other kind of problem, re-raise the exception
        else:
            raise e
    # Get the content from the object
    id_list_str = id_list_obj['Body'].read().decode('utf8').strip()
    id_str_list = id_list_str.splitlines()[args.start_index:args.end_index]
    random.shuffle(id_str_list)
    id_dict = get_id_dict([line.strip() for line in id_str_list])

    # Some combinations of options don't make sense:
    forbidden_combos = [('all', 'unread'), ('none', 'unread'), ('none', 'none')]
    assert (args.read_mode, args.stmt_mode) not in forbidden_combos, \
        ("The combination of reading mode %s and statement mode %s is not "
         "allowed." % (args.reading_mode, args.stmt_mode))

    # Init some timing dicts
    starts = {}
    ends = {}

    # Get a handle for the database
    if args.test:
        from indra.db import util as dbu
        db = dbu.get_test_db()
    else:
        db = None

    s3_log_prefix = ('reading_results/%s/logs/run_db_reading_queue/%s/'
                     % (args.basename, args.job_name))

    # Read everything ========================================
    starts['reading'] = datetime.now()
    outputs = produce_readings(id_dict, readers, verbose=True,
                               read_mode=args.read_mode,
                               get_preexisting=(args.stmt_mode == 'all'),
                               force_fulltext=args.force_fulltext,
                               prioritize=args.use_best_fulltext, db=db)
    ends['reading'] = datetime.now()

    # Preserve the sparser logs
    contents = os.listdir('.')
    sparser_logs = [fname for fname in contents
                    if fname.startswith('sparser') and fname.endswith('log')]
    sparser_log_dir = s3_log_prefix + 'sparser_logs/'
    for fname in sparser_logs:
        s3_key = sparser_log_dir + fname
        logger.info("Saving sparser logs to %s on s3 in %s."
                    % (s3_key, bucket_name))
        with open(fname, 'r') as f:
            client.put_object(Key=s3_key, Body=f.read(),
                              Bucket=bucket_name)

    # Convert the outputs to statements ==================================
    if args.stmt_mode != 'none':
        starts['statement production'] = datetime.now()
        stmt_data = produce_statements(outputs, n_proc=args.num_cores, db=db)
        ends['statement production'] = datetime.now()
    else:
        stmt_data = []

    report_statistics(outputs, stmt_data, starts, ends, s3_log_prefix, client,
                      bucket_name)

