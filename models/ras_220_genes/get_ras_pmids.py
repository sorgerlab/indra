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

# Load the lookup table of PMIDs with full texts (from the MySQL DB on EC2)
with open('pmids_fulltext.txt') as f:
    pmid_fulltexts = set([line.strip('\n') for line in f.readlines()])

with open('pmids_oa_txt.txt') as f:
    pmid_oa_txt = set([line.strip('\n') for line in f.readlines()])

with open('pmids_oa_xml.txt') as f:
    pmid_oa_xml = set([line.strip('\n') for line in f.readlines()])

with open('pmids_auth_xml.txt') as f:
    pmid_auth_xml = set([line.strip('\n') for line in f.readlines()])

pmid_map = {}
with open('pmid_pmcid_doi_map.txt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        doi = None if row[2] == 'None' else row[2]
        pmid_map[row[0]] = (row[1], doi)

def get_fulltexts(pmids_dict):
    fulltext_counts = OrderedDict()
    # Iterate over all 
    for gene, refs in pmids_dict.iteritems():
        refs_with_fulltext = pmid_fulltexts.intersection(set(refs))
        fulltext_counts[gene] = refs_with_fulltext
    return fulltext_counts

if __name__ == '__main__':
    import plot_formatting as pf

    (pmids, pmids_from_gene) = get_ids()

    not_found = 0
    ref_table = []
    for gene, refs in pmids_from_gene.iteritems():
        print gene
        for ref in refs:
            # Look up PMCID
            id_map = pmid_map.get(ref)
            if id_map is None:
                pmcid = None
                doi = None
            else:
                (pmcid, doi) = id_map
            # Look up full text status
            oa_xml = True if ref in pmid_oa_xml else False
            oa_txt = True if ref in pmid_oa_txt else False
            auth_xml = True if ref in pmid_auth_xml else False
            if pmcid is not None and not any((oa_xml, oa_txt, auth_xml)):
                not_found += 1
            row = (gene, ref, pmcid, doi, oa_xml, oa_txt, auth_xml)
            ref_table.append(row)

    import sys; sys.exit()


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
    pmids_gene_ft = get_fulltexts(pmids_from_gene)
    fig = plt.figure(figsize=(2, 2), dpi=300)
    ax = fig.gca()
    plot_parallel_counts(pmids_from_gene, pmids_gene_ft, ax)
    pf.format_axis(ax)
    plt.legend(['By Gene ID', 'Full Text in PMC'], loc='upper left',
               fontsize=pf.fontsize, frameon=False)
    plt.subplots_adjust(left=0.19, bottom=0.16)

    total_gene = np.sum([len(refs) for refs in pmids_from_gene.values()])
    total_gene_ft = np.sum([len(refs) for refs in pmids_gene_ft.values()])

    print "Total (by gene)", total_gene
    print "Total (by gene) with full text", total_gene_ft
    print "%.2f%% with full text" % \
          ((total_gene_ft / float(total_gene)) * 100)

    # Figure out how many of the publications have full text in PMC
    pmids_ft = get_fulltexts(pmids)
    fig = plt.figure(figsize=(2, 2), dpi=300)
    ax = fig.gca()
    plot_parallel_counts(pmids, pmids_ft, ax)
    pf.format_axis(ax)
    plt.legend(['By Name', 'Full Text in PMC'], loc='upper left',
               fontsize=pf.fontsize, frameon=False)
    plt.subplots_adjust(left=0.19, bottom=0.16)

    total_pmids = np.sum([len(refs) for refs in pmids.values()])
    total_pmids_ft = np.sum([len(refs) for refs in pmids_ft.values()])

    print "Total (by gene)", total_pmids
    print "Total (by gene) with full text", total_pmids_ft
    print "%.2f%% with full text" % \
          ((total_pmids_ft / float(total_pmids)) * 100)


    # Figure out how many of the publications have xml/txt/auth
    

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

