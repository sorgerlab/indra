import pickle
from indra.sources import trips

print('Loading text')
texts = pickle.load(open('text_subset.p', 'rb'))
print('\tDone!')

trips_statements = []
for text in texts:
    print(text)
    print('Parsing with TRIPS')
    tp = trips.process_text(text)
    print('\tDone!')
    trips_statements.extend(tp.statements)
    print(trips_statements)

