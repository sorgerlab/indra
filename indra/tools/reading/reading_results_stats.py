from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from matplotlib import pyplot as plt
import numpy as np
from indra.util import plot_formatting as pf
from collections import Counter
from indra.statements import *
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from copy import deepcopy
from indra.databases import uniprot_client
import sys
from indra.preassembler import grounding_mapper as gm

import logging

logger = logging.getLogger('reading_results')

pf.set_fig_params()


def load_file(stmts_file):
    logger.info("Loading results...")
    with open(stmts_file, 'rb') as f:
        results = pickle.load(f)
    return results


def report_stmt_counts(paper_stmts, plot_prefix=None):
    counts_per_paper = [(pmid, len(stmts)) for pmid, stmts in results.items()]
    counts = np.array([tup[1] for tup in counts_per_paper])
    logger.info("%.2f +/- %.2f statements per paper" % 
                (np.mean(counts), np.std(counts)))
    logger.info("Median %d statements per paper" % np.median(counts))

    zero_pmids = [pmid for pmid, stmts in results.items() if len(stmts) == 0]
    logger.info('%s papers with no statements' % len(zero_pmids))

    # Distribution of numbers of statements
    if plot_prefix:
        plot_filename = '%s_stmts_per_paper.pdf' % plot_prefix
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
        plt.savefig(plot_filename)


def plot_frequencies(counts, plot_filename, bin_interval):
    freq_dist = []
    bin_starts = xrange(0, len(counts), bin_interval)
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
        plot_filename = '%s_agent_distribution.pdf' % plot_prefix
        logger.info('Plotting occurrences of agent strings in %s' %
                    plot_filename)
        plot_frequencies(agent_counts, plot_filename, bin_interval)
    # Ungrounded agents
    ungrounded = gm.ungrounded_texts(stmts)
    logger.info('Top %d ungrounded strings, with frequencies' % list_length)
    for i in range(list_length):
        logger.info('%s: %d' % (ungrounded[i][0], ungrounded[i][1]))
    ungr_counts = [t[1] for t in ungrounded]
    if plot_prefix:
        plot_filename = '%s_ungrounded_distribution.pdf' % plot_prefix
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
    labels = [t.__name__ for t, n in sorted_counts]

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
    # Load the statements
    if len(sys.argv) < 2:
        print("Usage: %s reach_stmts_file" % sys.argv[0])
        sys.exit()
    results = load_file(sys.argv[1])
    all_stmts = [stmt for paper_stmts in results.values()
                      for stmt in paper_stmts]

    report_stmt_counts(results, plot_prefix='raw')
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

