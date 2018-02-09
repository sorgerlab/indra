from gensim.models.doc2vec import Doc2Vec
from gensim.models.doc2vec import TaggedDocument
from nltk.tokenize import sent_tokenize, word_tokenize

fname = 'abstracts_11000.txt'

# Read in sentences from file
sentences = []
with open(fname, 'r') as f:
    for line in f:
        sentences.append(line)

# Convert each sentence into a list of words
word_bags = [word_tokenize(s) for s in sentences]

alpha_word_bags = []
for i in range(len(word_bags)):
    awb = []
    for word in word_bags[i]:
        if word.isalnum():
            awb.append(word)
    alpha_word_bags.append(awb)

word_bags = alpha_word_bags

# Convert each list of words into a tagged document
tagged_docs = []
vocab = []
for i in range(len(word_bags)):
    vocab.extend(word_bags[i])
    tagged_docs.append( TaggedDocument(word_bags[i], (i,)) )

print(word_bags[0])
print(tagged_docs[0])

print(word_bags[1])
print(tagged_docs[1])

print('Now training the model')

# Train model
model = Doc2Vec(tagged_docs, size=50, min_count=2, workers=4, iter=100)


model.train(tagged_docs, total_examples=len(tagged_docs), epochs=model.iter)
model.save('paragraphVectors_abstract11000')


