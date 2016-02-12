import sys
from rdflib import Graph, Namespace, Literal
import csv

if __name__ == '__main__':
    rn = Namespace('indra/relations/')
    mn = Namespace('indra/modifications/')
    g = Graph()

    isa = rn.term('isa')

    g.add((mn.term('Phosphorylation'), isa, mn.term('Modification')))
    g.add((mn.term('PhosphorylationSerine'), isa, mn.term('Phosphorylation')))
    g.add((mn.term('PhosphorylationTyrosine'), isa, mn.term('Phosphorylation')))
    g.add((mn.term('PhosphorylationThreonine'), isa, mn.term('Phosphorylation')))

    with open('modification_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))
    
