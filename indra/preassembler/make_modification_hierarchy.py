import sys
from rdflib import Graph, Namespace, Literal
import csv

if __name__ == '__main__':
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    rn = Namespace(indra_ns + 'relations/')
    en = Namespace(indra_ns + 'entitiess/')
    g = Graph()

    isa = rn.term('isa')

    g.add((en.term('Phosphorylation'), isa, en.term('Modification')))
    g.add((en.term('PhosphorylationSerine'), isa, en.term('Phosphorylation')))
    g.add((en.term('PhosphorylationTyrosine'), isa, en.term('Phosphorylation')))
    g.add((en.term('PhosphorylationThreonine'), isa, en.term('Phosphorylation')))

    with open('modification_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))
    
