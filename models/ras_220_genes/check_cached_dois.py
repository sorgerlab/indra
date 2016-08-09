import csv
import pickle
from collections import Counter
import plot_formatting as pf
from matplotlib import pyplot as plt
from texttable import  Texttable

pf.set_fig_params()

pmid_map = {}
with open('pmid_pmcid_doi_map.txt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        doi = None if row[2] == '' else row[2]
        pmid_map[row[0]] = (row[1], doi)


print "Loading PMIDs from gene query"
with open('pmids_from_gene.pkl') as f:
    pmids = pickle.load(f)

doi_lookup = {}
print "Loading PMID->DOI cache"
with open('doi_cache.txt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        if not row[1]:
            continue
        doi_lookup[row[0]] = row[1]

print "Loading Pubmed metadata"
with open('pmids_from_gene_metadata.pkl') as f:
    pubmed_metadata = pickle.load(f)

total_pmids = [pmid for gene, pmid_list in pmids.iteritems()
                    for pmid in pmid_list]
unique_pmids = set(total_pmids)

pmids_no_doi = [pmid for pmid in unique_pmids \
                        if not (pmid_map.get(pmid) is not None and \
                                pmid_map.get(pmid)[1] is not None)]

print "Loading xref metadata"
with open('xref_metadata.pkl') as f:
    xref_metadata = pickle.load(f)

# Iterate over the PMIDs for which we didn't get DOIs from PubMed.
# Look each one up in the xref_metadata.
total = 0

not_in_doi_cache = []
no_xref_meta = []
no_pubmed_meta = []
matched = []
no_issn = []
no_page = []
page_mismatch = []
issn_mismatch = []

for pmid in pmids_no_doi:
    total += 1
    xref_doi = doi_lookup.get(pmid)
    if xref_doi is None:
        print pmid, "not found in doi_cache.txt"
        not_in_doi_cache.append(pmid)
        continue

    xref_meta = xref_metadata.get(xref_doi)
    if xref_meta is None:
        print pmid, xref_doi, "not found in xref_metadata.pkl"
        no_xref_meta.append(xref_doi)
        continue

    pubmed_meta = pubmed_metadata.get(pmid)
    if pubmed_meta is None:
        print pmid, "not found in pubmed_metadata"
        no_pubmed_meta.append(pmid)
        continue

    # ISSNs
    xr_issn_list = xref_meta.get('ISSN')
    pm_issn_list = pubmed_meta.get('issn_list')
    if not (xr_issn_list and pm_issn_list):
        print pmid, xref_doi, pm_issn_list, xr_issn_list, "Missing ISSNs"
        no_issn.append((pmid, xref_doi))
        continue
    matching_issns = set(pm_issn_list).intersection(set(xr_issn_list))

    # Page
    xref_page = xref_meta.get('page')
    pm_page = pubmed_meta.get('page')
    if not (xref_page and pm_page):
        print pmid, xref_doi, pm_page, xref_page, "Missing pages"
        no_page.append((pmid, xref_doi))
        continue
    pm_start = pm_page.split('-')[0].upper()
    xr_start = xref_page.split('-')[0].upper()
    if xr_start.endswith('E'):
        xr_start = 'E' + xr_start[:-1]

    # Check for match!
    if matching_issns and pm_start == xr_start:
        matched.append((pmid, xref_doi))
    elif not matching_issns:
        print "---ISSN Mismatch-----------"
        print pmid, xref_doi
        print "Pubmed ISSN:       ", pm_issn_list
        print "CrossR ISSN:       ", xr_issn_list
        print "Titles:"
        print "Pubmed", pubmed_meta.get('title')
        print "Xref", xref_meta.get('title')
        print pm_start, xr_start
        issn_mismatch.append((pmid, xref_doi))
    elif pm_start != xr_start:
        print "---Page Mismatch-----------"
        print pmid, xref_doi
        print "Pubmed ISSN:       ", pm_issn_list
        print "CrossR ISSN:       ", xr_issn_list
        print "Pages:", pm_start, xr_start
        print "Titles:"
        print "Pubmed", pubmed_meta.get('title')
        print "Xref", xref_meta.get('title')
        page_mismatch.append((pmid, xref_doi))
    else:
        assert False, "Should never happen"

print
print "Not in DOI cache:", len(not_in_doi_cache)
print "No XREF metadata:", len(no_xref_meta)
print "No Pubmed metadata:", len(no_pubmed_meta)
print "No ISSN:", len(no_issn)
print "No page:", len(no_page)
print "ISSN mismatch:", len(issn_mismatch)
print "Page mismatch:", len(page_mismatch)
print "Matched:", len(matched)

print "Collecting publishers"
publishers = []
for (pmid, doi) in matched:
    xref_meta = xref_metadata.get(doi)
    publishers.append((pmid, doi, xref_meta.get('publisher')))

has_license = []
for (pmid, doi) in matched:
    xref_meta = xref_metadata.get(doi)
    if xref_meta.get('license'):
        has_license.append((pmid, doi, xref_meta))

has_link = []
for (pmid, doi) in matched:
    xref_meta = xref_metadata.get(doi)
    if xref_meta.get('link'):
        has_link.append((pmid, doi, xref_meta))

jbc = []
for (pmid, doi) in matched:
    xref_meta = xref_metadata.get(doi)
    if xref_meta['publisher'] == 'American Society for Biochemistry & Molecular Biology (ASBMB)':
        jbc.append((pmid, doi, xref_meta))


pubs = [t[2] for t in publishers]
pubs_count = Counter(pubs)
pubs_count = sorted(pubs_count.items(), key=lambda x: x[1], reverse=True)
plt.ion()
fig = plt.figure(figsize=(2, 2), dpi=300)
ax = fig.gca()
ax.plot([t[1] for t in pubs_count])
ax.set_ylabel('No. of papers')
ax.set_xlabel('Publisher rank')
ax.set_xticks(range(0, 181, 30))
pf.format_axis(ax)
plt.subplots_adjust(left=0.22, bottom=0.16)
plt.show()
plt.savefig('publisher_distribution_linear.pdf')
ax.set_yscale('log')
plt.show()
plt.savefig('publisher_distribution_log.pdf')

# Make table of top 20
top_20_table = Texttable()
rows = [['Rank', 'Publisher']]
for i in range(20):
    rank = i + 1
    entry = '%s (%s)' % (str(pubs_count[i][0]), pubs_count[i][1])
    rows.append([rank, entry])
top_20_table.add_rows(rows)
print top_20_table.draw() + '\n'

