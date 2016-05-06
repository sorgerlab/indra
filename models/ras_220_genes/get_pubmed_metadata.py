from indra.literature import pubmed_client
import pickle

with open('pmids_from_gene.pkl') as f:
    pmids = pickle.load(f)

total_pmids = [pmid for gene, pmid_list in pmids.iteritems()
                    for pmid in pmid_list]
unique_pmids = sorted(list(set(total_pmids)), key=lambda x: int(x))

num_ids = len(unique_pmids)
chunk_size = 200
start_indices = range(0, num_ids, chunk_size)

results = {}
for start_ix in start_indices:
    print start_ix
    if start_ix + chunk_size < num_ids:
        end_ix = start_ix + chunk_size
    else:
        end_ix = num_ids
    results.update(pubmed_client.get_metadata_for_ids(
                                        unique_pmids[start_ix:end_ix]))
