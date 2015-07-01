import sys
import trips_client
from processor import TripsProcessor


def process_text(text):
    html = trips_client.send_query(text)
    xml = trips_client.get_xml(html)
    tp = TripsProcessor(xml)
    tp.get_complexes()
    tp.get_phosphorylation()
    tp.get_activating_mods()
    return tp

if __name__ == '__main__':
    input_fname = 'phosphorylate.xml'
    if len(sys.argv) > 1:
        input_fname = sys.argv[1]
    try:
        fh = open(input_fname, 'rt')
    except IOError:
        print 'Could not open file %s' % input_fname
        sys.exit()
    xml_string = fh.read()
    tp = TripsProcessor(xml_string)
    tp.get_complexes()
    tp.get_phosphorylation()
