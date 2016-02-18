import sys
from rdflib import Graph, Namespace, Literal
import csv

if __name__ == '__main__':
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    rn = Namespace(indra_ns + 'relations/')
    en = Namespace(indra_ns + 'entitiess/')
    g = Graph()
    
    isa = rn.term('isa')
    has_name = rn.term('hasName')

    g.add((en.term('Modification'), has_name, Literal('Modification')))
    g.add((en.term('Phosphorylation'), has_name, Literal('Phosphorylation')))
    g.add((en.term('Ubiquitination'), has_name, Literal('Ubiquitination')))
    g.add((en.term('Sumoylation'), has_name, Literal('Sumoylation')))
    g.add((en.term('Acetylation'), has_name, Literal('Acetylation')))
    g.add((en.term('Hydroxylation'), has_name, Literal('Hydroxylation')))
    g.add((en.term('PhosphorylationSerine'), has_name, 
            Literal('PhosphorylationSerine')))
    g.add((en.term('PhosphorylationTyrosine'), has_name, 
            Literal('PhosphorylationTyrosine')))
    g.add((en.term('PhosphorylationThreonine'), has_name, 
            Literal('PhosphorylationThreonine')))
    
    g.add((en.term('Phosphorylation'), isa, en.term('Modification')))
    g.add((en.term('Ubiquitination'), isa, en.term('Modification')))
    g.add((en.term('Acetylation'), isa, en.term('Modification')))
    g.add((en.term('Hydroxylation'), isa, en.term('Modification')))
    g.add((en.term('Sumoylation'), isa, en.term('Modification')))
    g.add((en.term('Phosphorylation'), isa, en.term('Modification')))
    g.add((en.term('PhosphorylationSerine'), isa, en.term('Phosphorylation')))
    g.add((en.term('PhosphorylationTyrosine'), isa, en.term('Phosphorylation')))
    g.add((en.term('PhosphorylationThreonine'), isa, en.term('Phosphorylation')))

    with open('modification_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))
    
