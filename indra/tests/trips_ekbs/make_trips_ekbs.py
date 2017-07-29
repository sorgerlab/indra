import re
from indra.processors import trips

def make_ekb(txt):
    fname = re.sub('[^a-zA-Z0-9]', '_', txt[:-1]) + '.ekb'
    tp = trips.process_text(txt, save_xml_name=fname)
    print(txt)
    print(tp.statements)

if __name__ == '__main__':
    lines = open('trips_ekb_sentences.txt', 'rt').readlines()
    for i, l in enumerate(lines):
        txt = l.strip()
        make_ekb(txt)
        print('---')
