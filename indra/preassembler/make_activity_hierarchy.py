from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from rdflib import Graph, Namespace, Literal
from os.path import abspath, dirname, join

hierarchy_path = join(dirname(abspath(__file__)),
                      '../resources/activity_hierarchy.rdf')

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

def main():
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    rn = Namespace(indra_ns + 'relations/')
    en = Namespace(indra_ns + 'entities/')
    g = Graph()

    isa = rn.term('isa')

    g.add((en.term('transcription'), isa, en.term('activity')))
    g.add((en.term('catalytic'), isa, en.term('activity')))
    g.add((en.term('gtpbound'), isa, en.term('activity')))
    g.add((en.term('kinase'), isa, en.term('catalytic')))
    g.add((en.term('phosphatase'), isa, en.term('catalytic')))
    g.add((en.term('gef'), isa, en.term('catalytic')))
    g.add((en.term('gap'), isa, en.term('catalytic')))

    save_hierarchy(g, hierarchy_path)

if __name__ == '__main__':
    main()
