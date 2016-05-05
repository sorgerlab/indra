import csv
from indra.literature import pubmed_client, crossref_client

doi_cache = {}
with open('doi_cache.txt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        doi_cache[row[0]] = row[1]

with open('missing_dois.txt') as f:
    missing_dois = [line.strip('\n') for line in f.readlines()]

def save(doi_cache, counter):
    with open('doi_cache_%.5d.txt' % counter, 'w') as f:
        print "Writing to doi cache"
        csvwriter = csv.writer(f, delimiter='\t')
        for k, v in doi_cache.iteritems():
            csvwriter.writerow((k, v))

for counter, ref in enumerate(missing_dois):
    if doi_cache.get(ref):
        continue
    title = pubmed_client.get_title(ref)
    if not title:
        print "No title, skipping", ref
        continue
    doi = crossref_client.doi_query(title)
    if doi:
        doi_cache[ref] = doi
        print "%d: %s --> %s" % (counter, ref, doi)
    else:
        print "No DOI for %s: %s" % (ref, title)
        continue

    if counter % 100 == 0:
        save(doi_cache, counter)

save(doi_cache, counter)
