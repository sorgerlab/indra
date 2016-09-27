from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
from rdflib import Graph, Namespace, Literal
import csv

if __name__ == '__main__':
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

    with open('../resources/activity_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))
