import sys
import trips_client
from processor import TripsProcessor


def process_text(text, save_xml_name='trips_output.xml'):
    html = trips_client.send_query(text)
    xml = trips_client.get_xml(html)
    if save_xml_name:
        trips_client.save_xml(xml, save_xml_name)
    return process_xml(xml)


def process_xml(xml_string):
    tp = TripsProcessor(xml_string)
    tp.get_complexes()
    tp.get_phosphorylation()
    tp.get_activating_mods()
    tp.get_activations()
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
