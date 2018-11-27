import sys
import json
from os.path import join, dirname, abspath
from rdflib import Graph, Namespace, Literal
from indra.sources import sofia


# Note that this is just a placeholder, it doesn't resolve as a URL
sofia_ns = Namespace('http://cs.cmu.edu/sofia/')
indra_ns = 'http://sorger.med.harvard.edu/indra/'
indra_rel_ns = Namespace(indra_ns + 'relations/')
isa = indra_rel_ns.term('isa')


def save_ontology(g, path):
    with open(path, 'wb') as out_file:
        g_bytes = g.serialize(format='nt')
        # Replace extra new lines in string and get rid of empty line at end
        g_bytes = g_bytes.replace(b'\n\n', b'\n').strip()
        # Split into rows and sort
        rows = g_bytes.split(b'\n')
        rows.sort()
        g_bytes = b'\n'.join(rows)
        out_file.write(g_bytes)


def build_ontology(ont_json, rdf_path):
    G = Graph()
    for top_key, entries in ont_json.items():
        for entry_key, examples in entries.items():
            if '/' in entry_key:
                parent, child = entry_key.split('/', maxsplit=1)
                parent_term = sofia_ns.term(parent)
                child_term = sofia_ns.term(entry_key)
                rel = (child_term, isa, parent_term)
                G.add(rel)
    save_ontology(G, rdf_path)


if __name__ == '__main__':
    # Path to a SOFIA ontology JSON file
    sofia_ont_json_file = sys.argv[1]
    with open(sofia_ont_json_file, 'r') as fh:
        sofia_ont_json = json.load(fh)
    sofia_rdf_path = join(dirname(abspath(sofia.__file__)),
                          'sofia_ontology.rdf')
    G = build_ontology(sofia_ont_json, sofia_rdf_path)
