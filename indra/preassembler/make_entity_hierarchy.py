from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
from rdflib import Graph, Namespace, Literal
import csv
from indra.util import read_unicode_csv
from os.path import join, dirname, abspath

hierarchy_path = join(dirname(abspath(__file__)),
                      '../resources/entity_hierarchy.rdf')

relations_file = join(dirname(abspath(__file__)),
                      '../resources/famplex/relations.csv')

hgnc_ns = Namespace('http://identifiers.org/hgnc.symbol/')
up_ns = Namespace('http://identifiers.org/uniprot/')
famplex_ns = Namespace('http://identifiers.org/fplx/')
indra_rel_ns = Namespace('http://sorger.med.harvard.edu/indra/relations/')

ns_map = {'http://identifiers.org/fplx/': 'FPLX',
          'http://sorger.med.harvard.edu/indra/relations/': 'FPLX',
          'http://identifiers.org/hgnc.symbol/': 'HGNC',
          'http://identifiers.org/uniprot/': 'UP',}

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

def make_term(ns_name, id):
    if ns_name == 'HGNC':
        term = hgnc_ns.term(id)
    elif ns_name == 'UP':
        term = up_ns.term(id)
    elif ns_name == 'FPLX':
        term = famplex_ns.term(id)
    else:
        raise ValueError("Unknown namespace %s" % ns_name)
    return term


def main(relations_file):
    g = Graph()
    isa = indra_rel_ns.term('isa')
    partof = indra_rel_ns.term('partof')

    family_names = set([])
    csv_rows = read_unicode_csv(relations_file, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL,
                                lineterminator='\r\n')
    for line in csv_rows:
        ns1, id1, rel, ns2, id2 = line
        term1 = make_term(ns1, id1)
        term2 = make_term(ns2, id2)
        if rel in ('isa', 'partof'):
            rel_term = indra_rel_ns.term(rel)
        else:
            raise ValueError("Invalid relation %s" % rel)
        g.add((term1, rel_term, term2))

    save_hierarchy(g, hierarchy_path)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        relations_file = sys.argv[1]
    main(relations_file)

