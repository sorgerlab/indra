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


def plot_parallel_counts(refs1, refs2, ax, labels, **kwargs):
    """Plot the distribution of reference counts, sorting both by the number
    of references for refs1."""
    pmid_counts = []
    for gene, pubs in refs1.iteritems():
        pmid_counts.append((gene, len(pubs), len(refs2[gene])))
    # Sort by the number of refs for the first arg list
    pmid_counts = sorted(pmid_counts, key=lambda x: x[1])
    ax.plot([x[2] for x in pmid_counts], color='b', label=labels[1], **kwargs)
    # Plot the sorted one last so it appears on top
    ax.plot([x[1] for x in pmid_counts], color='r', label=labels[0], **kwargs)
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
        doi = None if row[2] == '' else row[2]
        pmid_map[row[0]] = (row[1], doi)

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

    sys.exit()

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


    """
    doi_cache = {}
    with open('doi_cache.txt') as f:
        csvreader = csv.reader(f, delimiter='\t')
        for row in csvreader:
            doi_cache[row[0]] = row[1]

    total = 0
    no_text_or_doi = set([])
    no_cached_doi = set([])
    ref_table = []
    counter = 0
    for gene, refs in pmids_from_gene.iteritems():
        print gene
        for ref in refs:
            total += 1
            # Look up PMCID
            id_map = pmid_map.get(ref)
            if id_map is None:
                pmcid = None
                pm_doi = None
            else:
                (pmcid, pm_doi) = id_map
            # Look up full text status
            oa_xml = True if ref in pmid_oa_xml else False
            oa_txt = True if ref in pmid_oa_txt else False
            auth_xml = True if ref in pmid_auth_xml else False
            cached_doi = doi_cache.get(ref)
            if pm_doi and cached_doi:
                assert pm_doi == cached_doi
                print "DOIs match"
                doi = pm_doi
            elif pm_doi and not cached_doi:
                doi = pm_doi
                print "No cached DOI for", ref
                no_cached_doi.add(ref)
            elif cached_doi and not pm_doi:
                doi = cached_doi
            # Don't have DOI from anywhere
            elif not pm_doi and not cached_doi:
                title = pubmed_client.get_title(ref)
                if title:
                    print counter, "no doi for", ref
                    no_text_or_doi.add(ref)
                continue
                #    doi = crossref_client.doi_query(title)
                #    doi_cache[ref] = doi
                #    print "%d: Looked %s:%s --> %s" % (counter, gene, ref, doi)
                #    print title
            else:
                assert False #?????

            assert doi
            row = (gene, ref, pmcid, doi, oa_xml, oa_txt, auth_xml)
            ref_table.append(row)
            counter += 1

    print "Saving list of non-cached DOIs"
    with open('no_cached_doi.txt', 'w') as f:
        for ref in set(no_cached_doi):
            f.write('%s\n' % ref)

    # Remove duplicates by converting to a set
    with open('missing_dois.txt', 'w') as f:
        for ref in set(no_text_or_doi):
            f.write('%s\n' % ref)

    # Load whatever metadata we've got
    if os.path.isfile('xref_metadata.pkl'):
        with open('xref_metadata.pkl') as f:
            xref_meta = pickle.load(f)
    else:
        xref_meta = {}

    for counter, row in enumerate(ref_table):
        doi = row[3]
        # Do we already have metadata for this doi?
        if xref_meta.get(doi):
            print "Already have metadata for ", doi
            continue
        else:
            print "%d: querying for %s" % (counter, doi)
            metadata = crossref_client.get_metadata(doi)
            if metadata:
                xref_meta[doi] = metadata
            else:
                print "No metadata found for", doi
                continue
        if counter % 500 == 0:
            print "Saving metadata cache"
            with open('xref_metadata_%.5d.pkl' % counter, 'w') as f:
                pickle.dump(xref_meta, f)

    print "Final save of metadata cache"
    with open('xref_metadata_%.5d.pkl' % counter, 'w') as f:
        pickle.dump(xref_meta, f)

    import sys; sys.exit()

    # Randomly sample the PMIDs with no DOI to see if I can get the DOI
    # from PubMed
    #row_indices = range(len(ref_table))
    #sample_indices = np.random.choice(row_indices, size=100, replace=False)

    xr_found = []
    xr_not_found = []
    doi_cache = []
    import time
    start = time.time()
    for sample in samples:
        print "Querying for ", sample
        title = pubmed_client.get_title(sample)
        doi = crossref_client.doi_query(title)
        if doi:
            xr_found.append(sample)
            doi_cache.append((sample, doi))
        else:
            xr_not_found.append(sample)
    end = time.time()
    elapsed = end - start
    print "Elapsed time", elapsed

    with open('doi_cache.txt', 'w') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        csv.writerows(doi_cache)

    import sys; sys.exit()

    print "Querying for DOIs"
    pm_found = []
    pm_not_found = []
    for sample in samples:
        ids = id_lookup('PMID%s' % sample)
        if ids:
            if ids.get('doi'):
                pm_found.append(sample)
            else:
                pm_not_found.append(sample)
        else:
            pm_not_found.append(sample)
    """


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

