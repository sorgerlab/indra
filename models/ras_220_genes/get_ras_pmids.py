import csv
import sys
from indra.literature import pubmed_client
from collections import OrderedDict
from matplotlib import pyplot as plt
import pickle
import numpy as np
import boto3
import botocore
import os
from indra.literature import id_lookup
from indra.literature import crossref_client
from texttable import Texttable

def get_ids():
    """Search PubMed for references for the Ras 227 gene set."""
    # Check if we've got the files already
    if os.path.isfile('reading/pmids.pkl') and \
       os.path.isfile('reading/pmids_from_gene.pkl'):
        with open('reading/pmids.pkl') as pmids_file:
            pmids = pickle.load(pmids_file)
        with open('reading/pmids_from_gene.pkl') as pmids_from_gene_file:
            pmids_from_gene = pickle.load(pmids_from_gene_file)
        return (pmids, pmids_from_gene)

    # STEP 0: Get gene list
    gene_list = []
    # Get gene list from ras_pathway_proteins.csv
    with open('../../data/ras_pathway_proteins.csv') as f:
        csvreader = csv.reader(f, delimiter='\t')
        for row in csvreader:
            gene_list.append(row[0].strip())

    pmids = OrderedDict()
    pmids_from_gene = OrderedDict()

    for gene in gene_list:
        print "Querying for", gene
        ids_gene = set(pubmed_client.get_ids_for_gene(gene))
        print "Found %d in gene query" % len(ids_gene)
        # Hack to deal with excessive number of names
        if gene == 'MET':
            query_gene = 'CMET'
        elif gene == 'JUN':
            query_gene = 'CJUN'
        else:
            query_gene = gene
        ids_pubmed = set(pubmed_client.get_ids(query_gene,
                                               **{'retmax': 100000}))
        print "Found %d in string query" % len(ids_pubmed)
        pmids[gene] = ids_pubmed
        pmids_from_gene[gene] = ids_gene

    with open('reading/pmids.pkl', 'w') as f:
        pickle.dump(pmids, f)
    with open('reading/pmids_from_gene.pkl', 'w') as f:
        pickle.dump(pmids_from_gene, f)
    return (pmids, pmids_from_gene)


def plot_counts(refs, ax, **kwargs):
    """Plot the distribution of reference counts."""
    pmid_counts = []
    for gene, pubs in refs.iteritems():
        pmid_counts.append(len(pubs))
    pmid_counts = sorted(zip(refs.keys(), pmid_counts), key=lambda x: x[1])

    ax.plot(zip(*pmid_counts)[1], **kwargs)
    ax.set_yscale('log')
    ax.set_ylabel('Publications')
    ax.set_xlabel('Gene index')


def plot_parallel_counts(refs1, refs2, ax, labels, **kwargs):
    """Plot the distribution of reference counts, sorting both by the number
    of references for refs1."""
    pmid_counts = []
    for gene, pubs in refs1.iteritems():
        pmid_counts.append((gene, len(pubs), len(refs2[gene])))
    # Sort by the number of refs for the first arg list
    pmid_counts = sorted(pmid_counts, key=lambda x: x[1])
    ax.plot([x[1] for x in pmid_counts], color='r', label=labels[0],
            zorder=2, **kwargs)
    ax.plot([x[2] for x in pmid_counts], color='b', label=labels[1],
            zorder=1, **kwargs)
    ax.set_yscale('log')
    ax.set_ylabel('Publications')
    ax.set_xlabel('Gene index')
    return pmid_counts

"""
pmid_map = {}
with open('pmid_pmcid_doi_map.txt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        doi = None if row[2] == '' else row[2]
        pmid_map[row[0]] = (row[1], doi)
"""

def get_fulltexts(pmids_dict):
    fulltext_counts = OrderedDict()
    # Iterate over all 
    for gene, refs in pmids_dict.iteritems():
        refs_with_fulltext = pmid_fulltexts.intersection(set(refs))
        fulltext_counts[gene] = refs_with_fulltext
    return fulltext_counts

def num_unique_refs(pmids_dict):
    unique_refs = set([ref for gene, refs in pmids_dict.items()
                           for ref in refs])
    return len(unique_refs)

if __name__ == '__main__':
    import plot_formatting as pf

    # Load the PMIDs
    (pmids_by_name, pmids_by_gene) = get_ids()
    plt.ion()
    pf.set_fig_params()

    # Print some stats on the corpora
    # - Total and unique references returned
    total_name = np.sum([len(refs) for refs in pmids_by_name.values()])
    total_gene = np.sum([len(refs) for refs in pmids_by_gene.values()])
    refs_table = Texttable()
    refs_table.add_rows([['', 'By gene name', 'By gene ID'],
                    ['Total refs', total_name, total_gene],
                    ['Unique refs', num_unique_refs(pmids_by_name),
                                   num_unique_refs(pmids_by_gene)]])
    print refs_table.draw() + '\n'
    # Sort both reference lists by number of citations
    pmids_by_name_sorted = sorted([(gene, len(refs))
                                   for gene, refs in pmids_by_name.items()],
                                   key=lambda x: x[1], reverse=True)
    pmids_by_gene_sorted = sorted([(gene, len(refs))
                                   for gene, refs in pmids_by_gene.items()],
                                   key=lambda x: x[1], reverse=True)
    # - Top 10 genes for both searches
    top_10_table = Texttable()
    rows = [['Rank', 'By gene name (refs)', 'By gene ID (refs)']]
    for i in range(10):
        rank = i + 1
        by_name_str = '%s (%s)' % (pmids_by_name_sorted[i][0],
                                   pmids_by_name_sorted[i][1])
        by_gene_str = '%s (%s)' % (pmids_by_gene_sorted[i][0],
                                   pmids_by_gene_sorted[i][1])

        rows.append([rank, by_name_str, by_gene_str])
    top_10_table.add_rows(rows)
    print top_10_table.draw() + '\n'
    # - Bottom 10 genes for both searches
    bottom_10_table = Texttable()
    rows = [['Rank', 'By gene name (refs)', 'By gene ID (refs)']]
    for i in reversed(range(1, 11)):
        rank = len(pmids_by_name_sorted) - i + 1
        by_name_str = '%s (%s)' % (pmids_by_name_sorted[-i][0],
                                   pmids_by_name_sorted[-i][1])
        by_gene_str = '%s (%s)' % (pmids_by_gene_sorted[-i][0],
                                   pmids_by_gene_sorted[-i][1])
        rows.append([rank, by_name_str, by_gene_str])
    bottom_10_table.add_rows(rows)
    print bottom_10_table.draw() + '\n'

    # Plot citation distribution sorted by name search
    fig = plt.figure(figsize=(2, 2), dpi=300)
    ax = fig.gca()
    sorted_counts = plot_parallel_counts(pmids_by_name, pmids_by_gene, ax,
                                         labels=['By gene name', 'By gene ID'])
    plt.legend(loc='upper left', fontsize=pf.fontsize, frameon=False)
    pf.format_axis(ax)
    plt.subplots_adjust(left=0.19, bottom=0.16)
    plt.savefig('citations_by_gene.png', dpi=150)
    plt.savefig('citations_by_gene.pdf')

    # Figure out how many of the publications have full text in PMC
    dict_labels = ['By gene name', 'By gene ID']
    file_labels = ['pmids_by_name_ft', 'pmids_by_gene_ft']
    for dict_ix, pmids_dict in enumerate((pmids_by_name, pmids_by_gene)):
        pmids_dict_ft = get_fulltexts(pmids_dict)
        fig = plt.figure(figsize=(2, 2), dpi=300)
        ax = fig.gca()
        plot_parallel_counts(pmids_dict, pmids_dict_ft, ax,
                             labels=[dict_labels[dict_ix], 'In PMC OA'])
        pf.format_axis(ax)
        plt.legend(loc='upper left', fontsize=pf.fontsize, frameon=False)
        plt.subplots_adjust(left=0.19, bottom=0.16)
        plt.savefig('%s_line.png' % file_labels[dict_ix], dpi=150)
        plt.savefig('%s_line.pdf' % file_labels[dict_ix])

        print "Unique refs (%s): %s" % \
                (dict_labels[dict_ix], num_unique_refs(pmids_dict))
        print "Unique refs (%s) in PMC-OA: %s" % \
                (dict_labels[dict_ix], num_unique_refs(pmids_dict_ft))
        print "%.2f%% with full text" % \
              ((num_unique_refs(pmids_dict_ft) /
                  float(num_unique_refs(pmids_dict))) * 100)
        print
        # Get distribution of fractions in PMC
        pmids_dict_ft_fracs = []
        for gene in pmids_dict.keys():
            pmids_dict_ft_fracs.append((gene, len(pmids_dict_ft[gene]) /
                                          float(len(pmids_dict[gene]))))
        fig = plt.figure(figsize=(2, 2), dpi=300)
        ax = fig.gca()
        fracs = [tup[1] for tup in pmids_dict_ft_fracs]
        ax.hist(fracs, bins=np.linspace(0, 0.4, 20))
        ax.set_xlabel('Pct. PMC OA articles in search results')
        ax.set_xticks(np.linspace(0, 0.4, 5))
        ax.set_ylabel('Num of gene searches')
        plt.subplots_adjust(left=0.17, bottom=0.16)
        pf.format_axis(ax)
        plt.savefig('%s_hist.png' % file_labels[dict_ix], dpi=150)
        plt.savefig('%s_hist.pdf' % file_labels[dict_ix])
        # Expected fraction articles in PMC OA
        print "Mean %% in PMC OA: %s" % (np.mean(fracs) * 100)
        print "Stdev of %% in PMC OA: %s" % (np.std(fracs) * 100)
        print

