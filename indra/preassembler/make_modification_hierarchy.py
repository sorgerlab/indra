import sys
from rdflib import Graph, Namespace, Literal
import csv

if __name__ == '__main__':
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    rn = Namespace(indra_ns + 'relations/')
    en = Namespace(indra_ns + 'entities/')
    g = Graph()

    isa = rn.term('isa')
    has_name = rn.term('hasName')

    g.add((en.term('Modification'), has_name, Literal('modification')))
    g.add((en.term('Phosphorylation'), has_name, Literal('phosphorylation')))
    g.add((en.term('Ubiquitination'), has_name, Literal('ubiquitination')))
    g.add((en.term('Sumoylation'), has_name, Literal('sumoylation')))
    g.add((en.term('Acetylation'), has_name, Literal('acetylation')))
    g.add((en.term('Hydroxylation'), has_name, Literal('hydroxylation')))

    g.add((en.term('Phosphorylation'), isa, en.term('Modification')))
    g.add((en.term('Ubiquitination'), isa, en.term('Modification')))
    g.add((en.term('Sumoylation'), isa, en.term('Modification')))
    g.add((en.term('Acetylation'), isa, en.term('Modification')))
    g.add((en.term('Hydroxylation'), isa, en.term('Modification')))

    with open('modification_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))
