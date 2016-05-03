import csv
from indra.literature import pubmed_client
from collections import OrderedDict
from matplotlib import pyplot as plt
import pickle
import numpy as np
import boto3
import botocore
import os

def get_ids():
    """Search PubMed for references for the Ras 227 gene set."""
    # Check if we've got the files already
    if os.path.isfile('pmids.pkl') and os.path.isfile('pmids_from_gene.pkl'):
        with open('pmids.pkl') as pmids_file:
            pmids = pickle.load(pmids_file)
        with open('pmids_from_gene.pkl') as pmids_from_gene_file:
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
        ids_gene = pubmed_client.get_ids_for_gene(gene)
        print "Found %d in gene query" % len(ids_gene)
        # Hack to deal with excessive number of names
        if gene == 'MET':
            query_gene = 'CMET'
        elif gene == 'JUN':
            query_gene = 'CJUN'
        else:
            query_gene = gene
        ids_pubmed = pubmed_client.get_ids(query_gene, **{'retmax': 100000})
        print "Found %d in string query" % len(ids_pubmed)
        pmids[gene] = ids_pubmed
        pmids_from_gene[gene] = ids_gene

    with open('pmids.pkl', 'w') as f:
        pickle.dump(pmids, f)
    with open('pmids_from_gene.pkl', 'w') as f:
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


def plot_parallel_counts(refs1, refs2, ax, **kwargs):
    """Plot the distribution of reference counts, sorting both by the number
    of references for refs1."""
    pmid_counts = []
    for gene, pubs in refs1.iteritems():
        pmid_counts.append((gene, len(pubs), len(refs2[gene])))
    pmid_counts = sorted(pmid_counts, key=lambda x: x[1])
    ax.plot([x[1] for x in pmid_counts], color='r', **kwargs)
    ax.plot([x[2] for x in pmid_counts], color='b', **kwargs)
    ax.set_yscale('log')
    ax.set_ylabel('Publications')
    ax.set_xlabel('Gene index')
    return pmid_counts


if __name__ == '__main__':
    import plot_formatting as pf

    (pmids, pmids_from_gene) = get_ids()

    plt.ion()
    pf.set_fig_params()

    # Plot citation distribution for both methods
    fig = plt.figure(figsize=(2, 2), dpi=300)
    ax = fig.gca()
    plot_counts(pmids, ax, color='blue')
    plot_counts(pmids_from_gene, ax, color='red')
    pf.format_axis(ax)
    plt.subplots_adjust(left=0.19, bottom=0.16)
    plt.legend(['By Name', 'By Gene ID'], loc='upper left',
               fontsize=pf.fontsize, frameon=False)

    # Plot citation distribution sorted by name search
    fig = plt.figure(figsize=(2, 2), dpi=300)
    ax = fig.gca()
    sorted_counts = plot_parallel_counts(pmids, pmids_from_gene, ax)
    plt.legend(['By Name', 'By Gene ID'], loc='upper left',
               fontsize=pf.fontsize, frameon=False)
    pf.format_axis(ax)
    plt.subplots_adjust(left=0.19, bottom=0.16)

    # Figure out how many of the publications have full text in PMC
    # TODO
    import sys; sys.exit()

    s3 = boto3.resource('s3')
    bucket = s3.Bucket('bigmech')
    #http://stackoverflow.com/questions/33842944/check-if-a-key-exists-in-a-bucket-in-s3-using-boto3
    def is_in_s3(key_name):
        objs = list(bucket.objects.filter(Prefix=key_name))
        if len(objs) > 0: # and objs[0].key == key_name:
            print key_name, "found!"
            return True
        else:
            print key_name, "not found!"
            return False

    # Check open access for EGFR
    egfr_refs = pmids_from_gene['EGFR']
    found_refs = []
    for ref in egfr_refs:
        key_name = 'papers/PMID%s' % ref
        if is_in_s3(key_name):
            found_refs.append(ref)


    """
    pmid_counts = sorted(pmid_counts, key=lambda x: x[1])
    pub_counts = np.array([len(pmids[gene]) for gene in pmids.keys()])
    plt.plot(sorted(pub_counts,reverse=True))
    ax = plt.gca()
    ax.set_yscale('log')
    #plt.bar(pmids.keys(), 

    with open('pmid_list.tsv', 'wb') as f:
        csvwriter = csv.writer(f, delimiter='\t')

    """

