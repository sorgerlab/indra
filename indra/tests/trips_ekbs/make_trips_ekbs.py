import re
from indra import trips

lines = open('trips_ekb_sentences.txt', 'rt').readlines()
for i, l in enumerate(lines):
    txt = l.strip()
    fname = re.sub('[^a-zA-Z0-9]', '_', txt[:-1]) + '.ekb'
    tp = trips.process_text(txt, save_xml_name=fname)
    print(txt)
    print(tp.statements)
    print('---')
