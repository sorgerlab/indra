from indra.literature import *

# Pick an example gene
gene = 'SOS2'

# Get a list of PMIDs for the gene
pmids = pubmed_client.get_ids_for_gene(gene)

# Get the PMIDs that have XML in PMC
pmids_oa_xml = pmc_client.filter_pmids(pmids, 'oa_xml')

# Write the results to a file
with open('%s_pmids.txt' % gene, 'w') as f:
    for pmid in pmids_oa_xml:
        f.write('%s\n' % pmid)


