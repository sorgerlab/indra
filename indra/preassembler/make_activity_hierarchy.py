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

    g.add((en.term('Activity'), has_name, Literal('activity')))
    g.add((en.term('Kinase'), has_name, Literal('kinase')))
    g.add((en.term('Phosphatase'), has_name, Literal('phosphatase')))
    g.add((en.term('Catalytic'), has_name, Literal('catalytic')))
    g.add((en.term('GtpBound'), has_name, Literal('gtpbound')))
    g.add((en.term('Transcriptional'), has_name, Literal('transcriptional')))

    g.add((en.term('Transcriptional'), isa, en.term('Activity')))
    g.add((en.term('Catalytic'), isa, en.term('Activity')))
    g.add((en.term('GtpBound'), isa, en.term('Activity')))
    g.add((en.term('Kinase'), isa, en.term('Catalytic')))
    g.add((en.term('Phosphatase'), isa, en.term('Catalytic')))

    with open('activity_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))
