import os
import sys
import xml.etree.ElementTree as ET
from rdflib import Graph, Namespace, Literal


trips_ns = Namespace('http://trips.ihmc.us/concepts/')
indra_ns = 'http://sorger.med.harvard.edu/indra/'
indra_rel_ns = Namespace(indra_ns + 'relations/')
isa = indra_rel_ns.term('isa')


def save_hierarchy(g, path):
    with open(path, 'wb') as out_file:
        g_bytes = g.serialize(format='nt')
        # Replace extra new lines in string and get rid of empty line at end
        g_bytes = g_bytes.replace(b'\n\n', b'\n').strip()
        # Split into rows and sort
        rows = g_bytes.split(b'\n')
        rows.sort()
        g_bytes = b'\n'.join(rows)
        out_file.write(g_bytes)


def make_hierarchy(tree):
    g = Graph()
    concepts = tree.findall('concept')
    for concept in concepts:
        name = concept.attrib['name'].replace('ont::', '')
        if name == 'root':
            continue
        term = trips_ns.term(name)
        relations = concept.find("relation[@label='inherit']")
        related_names = [rr.strip().replace('ont::', '') for rr
                        in relations.text.strip().split('\n')]
        for related_name in related_names:
            related_term = trips_ns.term(related_name)
            g.add((term, isa, related_term))
    return g


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python make_trips_ontology.py /path/to/trips-ont-dsl.xml')
        sys.exit()
    fname = sys.argv[1]
    tree = ET.parse(fname)
    g = make_hierarchy(tree)
    save_hierarchy(g, 'trips_ontology.rdf')
