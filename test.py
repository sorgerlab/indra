from indra.sources.geneways.devel.find_pmid_xml import *

m25 = get_cached_full_text_mentions()

#s = FullTextMentionSet(m25)
#s.collect_tag_names()
#print(s.tag_names)

strip_tags = ['italic', 'bold', 'sup', 'sub', 'xref']
block_tags = ['OAI-PMH', 'responseDate', 'request', 'header', 'journal-meta', 'back', 'fig', 'formula', 'id']
remove_tags = ['ext-link']

idx = 2

t = m25[0].xml_full_text
f = open('debug.txt', 'w')
f.write(t)
f.close()

for idx in range(len(m25)):
    matches = m25[idx].find_matching_sentences()
    print(str(idx) + ": " + repr(m25[idx]))
    for sentence in matches:
        print('\t* %s' % sentence)
    print('\n')
