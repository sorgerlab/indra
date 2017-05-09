"""
Script to get curated PMIDs in Entrez Gene for all approved gene symbols in
HGNC. There are roughly 40k approved symbols so the search takes several hours.
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

if __name__ == '__main__':
    import os
    import pickle
    from indra.util import read_unicode_csv
    from indra.literature import pubmed_client

    # Get all HGNC IDs
    hgnc_file = '../../indra/resources/hgnc_entries.tsv'
    lines = read_unicode_csv(hgnc_file, delimiter='\t')
    # Skip the header line
    next(lines)
    hgnc_names = [line[1] for line in lines
                          if line[3] == 'Approved']
    # Load the dict of PMIDs, if there is one
    dict_filename = 'pmids_for_gene.pkl'
    if os.path.exists(dict_filename):
        with open(dict_filename, 'rb') as f:
            pmids_for_gene = pickle.load(f)
    else:
        pmids_for_gene = {}
    # Get PMIDs for each HGNC ID
    num_added = 0
    for hgnc_name in hgnc_names:
        # If HGN
        #print('Getting PMIDs for %s' % hgnc_name)
        if hgnc_name in pmids_for_gene:
            print('%s: already got PMIDs, skipping' % hgnc_name)
            continue
        try:
            pmids = pubmed_client.get_ids_for_gene(hgnc_name)
        except ValueError as ex:
            print("Exception in gettting PMIDs for %s: %s" % (hgnc_name, ex))
            print("Continuing...")
            continue
        print('%s: %d PMIDs' % (hgnc_name, len(pmids)))
        pmids_for_gene[hgnc_name] = pmids
        num_added += 1
        if num_added % 50 == 0:
            print("Saving info for %d genes" % len(pmids_for_gene))
            with open(dict_filename, 'wb') as f:
                pickle.dump(pmids_for_gene, f)
            unique_pmids = set([pmid for pmid_list in pmids_for_gene.values()
                                     for pmid in pmid_list])
            print('Total PMIDs: %d' % len(unique_pmids))
