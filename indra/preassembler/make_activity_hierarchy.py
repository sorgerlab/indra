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
    act = Namespace(indra_ns + 'activities/')
    g = Graph()

    isa = rn.term('isa')

    g.add((act.term('transcription'), isa, act.term('activity')))
    g.add((act.term('catalytic'), isa, act.term('activity')))
    g.add((act.term('gtpbound'), isa, act.term('activity')))
    g.add((act.term('kinase'), isa, act.term('catalytic')))
    g.add((act.term('phosphatase'), isa, act.term('catalytic')))
    g.add((act.term('gef'), isa, act.term('catalytic')))
    g.add((act.term('gap'), isa, act.term('catalytic')))

    save_hierarchy(g, hierarchy_path)

if __name__ == '__main__':
    main()
