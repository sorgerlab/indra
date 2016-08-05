from indra.literature import *
from functools import partial

# Pick an example gene
gene = 'KRAS'

# Get a list of PMIDs for the gene
pmids = pubmed_client.get_ids_for_gene(gene)

# Get the PMIDs that have XML in PMC
pmids_oa_xml = pmc_client.filter_pmids(pmids, 'oa_xml')
pmids_auth_xml = pmc_client.filter_pmids(pmids, 'auth_xml')



"""
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))
print len(pmids)

# Get PMC/DOI and other metadata for the list of papers, 200 at a time
results = map(pubmed_client.get_metadata_for_ids, chunker(pmids, 200))
"""
