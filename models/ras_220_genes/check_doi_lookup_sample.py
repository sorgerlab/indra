import numpy as np
import csv
from indra.literature import pubmed_client, crossref_client

pmid_map = {}
with open('pmid_pmcid_doi_map.txt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        doi = None if row[2] == '' else row[2]
        pmid_map[row[0]] = (row[1], doi)

with open('no_cached_doi.pkl') as f:
    no_cached_doi = [line.strip('\n') for line in f.readlines()]

# Get random sample of 100 non-cached DOIs
row_indices = range(len(no_cached_doi))
sample_indices = np.random.choice(row_indices, size=100, replace=False)

total_tests = 0
num_passes = 0
for sample_ix in sample_indices:
    ref = no_cached_doi[sample_ix]
    pm_doi = pmid_map[ref][1]
    xref_doi = crossref_client.doi_query(ref)
    if xref_doi:
        total_tests += 1
        if xref_doi.lower() == pm_doi.lower():
            num_passes += 1
            print "Pass: %s / %s" % (num_passes, total_tests)
        else:
            print "Fail: pm: %s, xref: %s" % (pm_doi, xref_doi)
    else:
        print "No DOI found, skipping"

