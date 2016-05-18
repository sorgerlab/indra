import csv
import pickle
from indra.literature import pubmed_client

pmid_map = {}
with open('pmid_pmcid_doi_map.txt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        doi = None if row[2] == '' else row[2]
        pmid_map[row[0]] = (row[1], doi)

print "Loading xref metadata"
with open('xref_metadata.pkl') as f:
    xref_metadata = pickle.load(f)

#print "Loading Pubmed metadata"
#with open('pmids_from_gene_metadata.pkl') as f:
#    pubmed_metadata = pickle.load(f)

print "Loading PMIDs from gene query"
with open('pmids_from_gene.pkl') as f:
    pmids = pickle.load(f)

total_pmids = [pmid for gene, pmid_list in pmids.iteritems()
                    for pmid in pmid_list]
unique_pmids = sorted(list(set(total_pmids)), key=lambda x: int(x))

# Iterate over the PMIDs in the list
counter = 0
mismatch = []
for pmid in unique_pmids:
    counter += 1
    pmid_result = pmid_map.get(pmid)
    if pmid_result and pmid_result[1]:
        doi = pmid_result[1]
        # Lookup the metadata for the pub in pubmed
        #pubmed_meta = pubmed_metadata[pmid]
        pubmed_meta = pubmed_client.get_metadata_for_ids([pmid])[pmid]
        # Look up the metadata for the pub in CrossRef
        xref_meta = xref_metadata.get(doi)
        if xref_meta is None:
            print counter, pmid, doi, "Not found in Xref, skipping"
            continue
        xr_issn_list = xref_meta.get('ISSN')
        if xr_issn_list is None:
            print counter, pmid, doi, "No ISSNs in XREF, skipping"
            continue
        # PM ISSN
        nlm_id = pubmed_meta['journal_nlm_id']
        pm_issn = pubmed_meta['issn']
        pm_issn_linking = pubmed_meta['issn_linking']
        pm_issn_list = pubmed_client.get_issns_for_journal(nlm_id)
        if pm_issn_list is None:
            pm_issn_list = []
        pm_issn_list += [pm_issn, pm_issn_linking]

        #print counter, pmid, nlm_id, pubmed_meta['issn'], xr_issn_list, \
        #      "No ISSNs in PM, skipping"
        #continue

        matching_issns = set(pm_issn_list).intersection(set(xr_issn_list))
        if not matching_issns:
            print counter, "---Mismatch-----------"
            print "Pubmed ISSN:       ", pm_issn_list
            print "CrossR ISSN:       ", xr_issn_list
            mismatch.append((pm_issn_list, xr_issn_list))
        else:
            xref_page = xref_meta.get('page')
            pm_page = pubmed_meta.get('page')
            if xref_page is None or pm_page is None:
                print counter, pmid, doi, "no page", pm_page, xref_page
                continue
            #pm_page_expanded = pubmed_client.expand_pagination(pm_page)
            pm_start = pm_page.split('-')[0].upper()
            xr_start = xref_page.split('-')[0].upper()
            if xr_start.endswith('E'):
                xr_start = 'E' + xr_start[:-1]
            if xr_start == pm_start:
                pass
            else:
                print counter, pmid, "issn match, page mismatch", \
                      pm_start, xr_start
            #print counter, "--- Match!!! --------"
            #print "Pubmed ISSN:", pm_issn_list
            #print "CrossR ISSN:", xr_issn_list
            #print "Match      :", matching_issns
    else:
        pass
        #print "No PM entry or doi for", pmid, "skipping"
