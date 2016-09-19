import sys
from rdflib import Graph, Namespace, Literal
import csv

if __name__ == '__main__':
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    rn = Namespace(indra_ns + 'relations/')
    en = Namespace(indra_ns + 'entities/')
    g = Graph()

    isa = rn.term('isa')

    g.add((en.term('phosphorylation'), isa, en.term('modification')))
    g.add((en.term('ubiquitination'), isa, en.term('modification')))
    g.add((en.term('sumoylation'), isa, en.term('modification')))
    g.add((en.term('acetylation'), isa, en.term('modification')))
    g.add((en.term('hydroxylation'), isa, en.term('modification')))

    with open('../resources/modification_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))
