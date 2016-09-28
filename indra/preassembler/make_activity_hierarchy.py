from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from rdflib import Graph, Namespace, Literal
from os.path import abspath, dirname, join

hierarchy_path = join(dirname(abspath(__file__)),
                      '../resources/activity_hierarchy.rdf')

def main():
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    rn = Namespace(indra_ns + 'relations/')
    en = Namespace(indra_ns + 'entities/')
    g = Graph()

    isa = rn.term('isa')

    g.add((en.term('transcriptional'), isa, en.term('activity')))
    g.add((en.term('catalytic'), isa, en.term('activity')))
    g.add((en.term('gtpbound'), isa, en.term('activity')))
    g.add((en.term('kinase'), isa, en.term('catalytic')))
    g.add((en.term('phosphatase'), isa, en.term('catalytic')))

    with open(hierarchy_path, 'wb') as out_file:
        g_bytes = g.serialize(format='xml', encoding='utf-8')
        out_file.write(g_bytes)

if __name__ == '__main__':
    main()
