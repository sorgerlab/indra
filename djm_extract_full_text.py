from indra.literature.pubmed_client import get_abstract
import pickle
import random
from nltk.tokenize import sent_tokenize, word_tokenize
import codecs

l = pickle.load(open('pmids_multiSource.p', 'rb'))
pmids = list(l['reach_statement_list.p'])
random.shuffle(pmids)

with codecs.open('abstracts_500.txt', 'w', encoding='utf-8') as f:
    for i in range(500):
        print(i)
        pmid = pmids[i]
        try:
            a = get_abstract(pmid, False)

            sentences = sent_tokenize(a)

            if i > 0:
                f.write('\n')
            f.write('\n'.join(sentences))
        except:
            pass

