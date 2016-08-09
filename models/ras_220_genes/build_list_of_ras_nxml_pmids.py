import pickle

print "Loading PMIDs from gene query"
with open('pmids_from_gene.pkl') as f:
        pmids = pickle.load(f)
total_pmids = [pmid for gene, pmid_list in pmids.iteritems()
                    for pmid in pmid_list]
unique_pmids = set(total_pmids)

# Load PMC ids with XML
with open('pmids_oa_xml.txt') as f:
    pmids_oa_xml = set([line.strip('\n') for line in f.readlines()]) 

with open('pmids_auth_xml.txt') as f:
    pmids_auth_xml = set([line.strip('\n') for line in f.readlines()]) 

pmids_with_oa = unique_pmids.intersection(pmids_oa_xml)
pmids_with_auth = unique_pmids.intersection(pmids_auth_xml)
pmids_any_xml = pmids_with_oa.union(pmids_with_auth)

with open('ras_gene_pmids_oa_xml.txt', 'w') as f:
    for pmid in pmids_with_oa:
        f.write('%s\n' % pmid)

with open('ras_gene_pmids_auth_xml.txt', 'w') as f:
    for pmid in pmids_with_auth:
        f.write('%s\n' % pmid)

with open('ras_gene_pmids_any_xml.txt', 'w') as f:
    for pmid in pmids_any_xml:
        f.write('%s\n' % pmid)

