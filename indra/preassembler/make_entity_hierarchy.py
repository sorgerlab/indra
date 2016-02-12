import sys
from rdflib import Graph, Namespace, Literal
import csv

if __name__ == '__main__':
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    if len(sys.argv) > 1:
        proteins_file = sys.argv[1]
    else:
        proteins_file = '../../data/ras_pathway_proteins.csv'
    rn = Namespace(indra_ns + 'relations/')
    en = Namespace(indra_ns + 'entities/')
    g = Graph()

    has_name = rn.term('hasName')
    has_synonym = rn.term('hasSynonym')
    isa = rn.term('isa')

    with open(proteins_file) as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            symbol, name, synonyms_str, family = line
            g.add((en.term(symbol), has_name, Literal(name)))
            synonyms = synonyms_str.split(', ')
            for s in synonyms:
                g.add((en.term(symbol), has_synonym, Literal(s)))
            if family.strip() != '':
                g.add((en.term(symbol), isa, en.term(family)))
    with open('entity_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))
    
