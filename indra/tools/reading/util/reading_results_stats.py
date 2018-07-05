from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import pickle
import logging
import numpy as np
from datetime import datetime
from collections import Counter
from matplotlib import pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

logger = logging.getLogger('reading_results')


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Get statistics on a bunch of statements.'
        )

    subparsers = parser.add_subparsers(title='Source')
    subparsers.required = True
    subparsers.dest = 'source'

    # Make file subparser
    file_subparser_parent = ArgumentParser(add_help=False)
    file_subparser_parent.add_argument(
        'file_path',
        help=('The path to the pickle file containing the statements to be '
              'analyzed.')
        )
    file_parser = subparsers.add_parser(
        'from-pickle',
        parents=[file_subparser_parent],
        description=('Get statistics of statements in a pickle file.'),
        formatter_class=ArgumentDefaultsHelpFormatter
        )

    # Make db subparser
    db_subparser_parent = ArgumentParser(add_help=False)
    db_subparser_parent.add_argument(
        '--indra_version',
        help='Specify the indra version for the batch of statements'
        )
    db_subparser_parent.add_argument(
        '--date_range',
        help=('Specify the range of datetimes for statements. Must be in the '
              'format: \"YYYYMMDDHHMMSS:YYYMMDDHHMMSS\". If you do not want '
              'to impose the upper or lower bound, simply leave it blank, eg. '
              '\"YYYYMMDDHHMMSS:\" if you don\'t care about the upper bound.')
        )
    db_parser = subparsers.add_parser(
        'from-db',
        parents=[db_subparser_parent],
        description=('Get statistics from statements on the database.'),
        formatter_class=ArgumentDefaultsHelpFormatter
        )

    # Parse the arguments.
    args = parser.parse_args()


from indra.util import plot_formatting as pf
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler import grounding_mapper as gm

#pf.set_fig_params()

# TODO: This file will definitely need some fixing.


def load_file(stmts_file):
    logger.info("Loading results...")
    with open(stmts_file, 'rb') as f:
        results = pickle.load(f)
    return results


def load_stmts_from_db(clauses, db):
    assert clauses, "No constraints provided for selecting statements."

    stmt_query = db.filter_query(
        [db.Statements, db.Readings.reader, db.TextRef],
        *(clauses + db.join(db.RawStatements, db.TextRef))
        )
    all_db_stmt_data = stmt_query.all()
    logger.info("Found %d statements on the database." % len(all_db_stmt_data))
    all_stmt_tr_pairs = [(_get_statement_object(db_stmt), tr)
                         for db_stmt, _, tr in all_db_stmt_data]
    map(_set_evidence_text_ref, all_stmt_tr_pairs)
    all_stmts = [stmt for stmt, _ in all_stmt_tr_pairs]

    results_db = {'reach': {}, 'sparser': {}}
    for db_stmt, reader_name, tr in all_db_stmt_data:
        reader_name = reader_name.lower()
        stmt = _get_statement_object(db_stmt)
        _set_evidence_text_ref(stmt, tr)
        if tr.id in results_db[reader_name]:
            results_db[reader_name][tr.id].append(stmt)
        else:
            results_db[reader_name][tr.id] = [stmt]

    results = {reader_name: {trid: stmts for trid, stmts in paper_stmts.items()}
               for reader_name, paper_stmts in results_db.items()}
    logger.info("Done sorting db statements.")
    return all_stmts, results


def report_stmt_counts(results, plot_prefix=None):
    counts_per_paper = [(pmid, len(stmts)) for pmid, stmts in results.items()]
    counts = np.array([tup[1] for tup in counts_per_paper])
    logger.info("%.2f +/- %.2f statements per paper" %
                (np.mean(counts), np.std(counts)))
    logger.info("Median %d statements per paper" % np.median(counts))

    zero_pmids = [pmid for pmid, stmts in results.items() if len(stmts) == 0]
    logger.info('%s papers with no statements' % len(zero_pmids))

    # Distribution of numbers of statements
    if plot_prefix:
        plot_filename = '%s_stmts_per_paper.png' % plot_prefix
        logger.info('Plotting distribution of statements per paper in %s' %
                    plot_filename)
        fig = plt.figure(figsize=(2, 2), dpi=300)
        ax = fig.gca()
        ax.hist(counts, bins=20)
        #ax.set_xticks([0, 100, 200, 300, 400])
        pf.format_axis(ax)
        plt.subplots_adjust(left=0.23, bottom=0.16)
        ax.set_xlabel('No. of statements')
        ax.set_ylabel('No. of papers')
        fig = plt.gcf()
        fig.savefig(plot_filename, format='png')
        logger.info('plotted...')


def plot_frequencies(counts, plot_filename, bin_interval):
    logger.info("Plotting frequencies")
    freq_dist = []
    bin_starts = range(0, len(counts), bin_interval)
    for bin_start_ix in bin_starts:
        bin_end_ix = bin_start_ix + bin_interval
        if bin_end_ix < len(counts):
            freq_dist.append(np.sum(counts[bin_start_ix:bin_end_ix]))
        else:
            freq_dist.append(np.sum(counts[bin_start_ix:]))
    freq_dist = np.array(freq_dist)
    fracs_total = np.cumsum(freq_dist) / float(np.sum(counts))

    fig = plt.figure(figsize=(2, 2), dpi=300)
    ax = fig.gca()
    ax.plot(bin_starts, fracs_total)
    pf.format_axis(ax)
    plt.subplots_adjust(left=0.23, bottom=0.16)
    ax.set_xlabel('String index')
    ax.set_ylabel('No. of occurrences')
    ax.set_ylim([0, 1])
    plt.savefig(plot_filename)


def report_grounding(stmts, list_length=10, bin_interval=10, plot_prefix=None):
    # All agents
    agents = gm.agent_texts_with_grounding(stmts)
    logger.info('Top %d agent strings, with frequencies:' % list_length)
    for i in range(list_length):
        logger.info('%s: %d' % (agents[i][0], agents[i][2]))
    agent_counts = [t[2] for t in agents]
    if plot_prefix:
        plot_filename = '%s_agent_distribution.png' % plot_prefix
        logger.info('Plotting occurrences of agent strings in %s' %
                    plot_filename)
        plot_frequencies(agent_counts, plot_filename, bin_interval)
    # Ungrounded agents
    ungrounded = gm.ungrounded_texts(stmts)
    logger.info('Top %d ungrounded strings, with frequencies' % list_length)
    for i in range(min(len(ungrounded), list_length)):
        logger.info('%s: %d' % (ungrounded[i][0], ungrounded[i][1]))
    ungr_counts = [t[1] for t in ungrounded]
    if plot_prefix:
        plot_filename = '%s_ungrounded_distribution.png' % plot_prefix
        logger.info('Plotting occurrences of ungrounded agents in %s' %
                    plot_filename)
        plot_frequencies(ungr_counts, plot_filename, bin_interval)

    # What fraction of statements grounded?
    all_ungrounded = 0
    any_ungrounded = 0
    for stmt in stmts:
        agents_ungrounded = []
        for ag in stmt.agent_list():
            if ag is not None and ag.db_refs.keys() == ['TEXT']:
                agents_ungrounded.append(True)
            else:
                agents_ungrounded.append(False)
        if all(agents_ungrounded):
            all_ungrounded += 1
        if any(agents_ungrounded):
            any_ungrounded += 1

    logger.info('%d of %d (%.1f%%) statements with all agents ungrounded' %
                (all_ungrounded, len(stmts),
                 100 * (all_ungrounded / float(len(stmts)))))
    logger.info('%d of %d (%.1f%%) statements with any agents ungrounded' %
                (any_ungrounded, len(stmts),
                 100 * (any_ungrounded / float(len(stmts)))))


def report_stmt_types(stmts, plot_prefix=None):
    # Number of statements of different types
    stmt_types = [type(stmt) for stmt in stmts]
    stmt_type_counter = Counter(stmt_types)

    fig = plt.figure(figsize=(2, 3), dpi=300)
    ax = fig.gca()
    sorted_counts = sorted(stmt_type_counter.items(), key=lambda x: x[1],
                           reverse=True)
    labels = [t.__name__ for t, _ in sorted_counts]

    logger.info('Distribution of statement types:')
    for stmt_type, count in sorted_counts:
        logger.info('%s: %d' % (stmt_type.__name__, count))

    if plot_prefix:
        x_coords = np.arange(len(sorted_counts))
        ax.bar(x_coords, [x[1] for x in sorted_counts])
        ax.set_xticks(x_coords + 0.4)
        ax.set_xticklabels(labels, rotation=90)
        pf.format_axis(ax)
        ax.set_ylabel('No. of statements')
        plt.subplots_adjust(left=0.29, bottom=0.34)
        plt.savefig('%s_stmt_types.pdf' % plot_prefix)


def report_stmt_participants(stmts):
    missing_participants = 0
    for stmt in stmts:
        if any([True if ag is None else False for ag in stmt.agent_list()]):
            missing_participants += 1
    logger.info('%d of %d (%.1f%%) statements missing participants' %
                (missing_participants, len(stmts),
                 100 * (missing_participants / float(len(stmts)))))


def report_evidence_distribution(stmts, list_length=10, plot_prefix=None):
    # Sort by evidence
    sorted_stmts = sorted(stmts, key=lambda x: len(x.evidence), reverse=True)
    logger.info('Top %d statements by evidence:' % list_length)
    for i in range(list_length):
        logger.info('%s: %d' % (sorted_stmts[i], len(sorted_stmts[i].evidence)))

    # Distribution of pieces of evidence
    if plot_prefix:
        fig = plt.figure(figsize=(2, 2), dpi=300)
        ax = fig.gca()
        ax.plot([len(stmt.evidence) for stmt in sorted_stmts])
        pf.format_axis(ax)
        ax.set_xlabel('Statement index')
        ax.set_ylabel('No. of sentences')
        ax.set_yscale('log')
        plt.subplots_adjust(left=0.23, bottom=0.16)
        plt.savefig('%s_evidence_dist.pdf' % plot_prefix)
    return sorted_stmts


if __name__ == '__main__':
    # Load the statements.
    if args.source == 'from-pickle':
        logger.info("Getting statements from pickle file.")
        results = load_file(args.file_path)
        all_stmts = [stmt for reader_stmts in results.values()
                     for paper_stmts in reader_stmts
                     for stmt in paper_stmts]
    if args.source == 'from-db':
        logger.info("Getting statements from the database.")
        from indra.db import get_primary_db
        from indra.db.util import get_raw_stmts_frm_db_list, \
            _get_statement_object, _set_evidence_text_ref

        db = get_primary_db()
        clauses = []
        if args.indra_version:
            clauses.append(db.RawStatements.indra_version == args.indra_version)
        if args.date_range:
            min_date_str, max_date_str = args.date_range.split(':')
            if min_date_str:
                min_date = datetime.strptime(min_date_str, '%Y%m%d%H%M%S')
                clauses.add(db.RawStatements.create_date > min_date)
            if max_date_str:
                max_date = datetime.strptime(max_date_str, '%Y%m%d%H%M%S')
                clauses.add(db.RawStatements.create_date < max_date)

        all_stmts, results = load_stmts_from_db(clauses, db)

    report_stmt_counts(results['reach'], plot_prefix='raw_reach')
    report_stmt_counts(results['sparser'], plot_prefix='raw_sparser')
    report_grounding(all_stmts, plot_prefix='raw')
    report_stmt_types(all_stmts, plot_prefix='raw')
    report_stmt_participants(all_stmts)

    # Map grounding
    logger.info('Mapping grounding...')
    gmap = gm.GroundingMapper(gm.default_grounding_map)
    map_stmts = gmap.map_agents(all_stmts)

    report_grounding(map_stmts, plot_prefix='preassembled')

    # Combine duplicates
    logger.info('Removing duplicates...')
    pa = Preassembler(hierarchies, map_stmts)
    pa.combine_duplicates()

    report_evidence_distribution(pa.unique_stmts, plot_prefix='preassembled')
