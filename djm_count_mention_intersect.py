from indra.sources.geneways.geneways_actionmention_parser import *

def get_pmids(parser):
    pmids = []
    mentions = parser.action_mentions
    for mention in mentions:
        pmids.append(mention.pmid)
    return set(pmids)

fname_old = '/Users/daniel/data/geneways_data/ivan_original_files/human_actionmention.txt'
fname_new = '/Users/daniel/data/geneways_data/human_db_dump_21jun15/human_actionmention.txt'

parser_old = GenewaysActionMentionParser(fname_old)
parser_new = GenewaysActionMentionParser(fname_new)

pmids_old = get_pmids(parser_old)
pmids_new = get_pmids(parser_new)

inter = pmids_old.intersection(pmids_new)
u = pmids_old.union(pmids_new)

print('Number of old PMIDs:', len(pmids_old))
print('Number of new PMIDs:', len(pmids_new))
print('Union size:', len(u))
print('Intersection size:', len(inter))
