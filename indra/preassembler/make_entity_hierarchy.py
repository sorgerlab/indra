import sys
from rdflib import Graph, Namespace, Literal
import csv
import urllib2


def make_term(ns_name, id):
    if ns_name == 'HGNC':
        term = hgnc_ns.term(id)
    elif ns_name == 'UP':
        term = up_ns.term(id)
    elif ns_name == 'BE':
        term = en.term(id)
    else:
        raise ValueError("Unknown namespace %s" % ns)
    return term


if __name__ == '__main__':
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    if len(sys.argv) > 1:
        relations_file = sys.argv[1]
    else:
        relations_file = '../../bioentities/relations.csv'
    rn = Namespace(indra_ns + 'relations/')
    en = Namespace(indra_ns + 'entities/')
    hgnc_ns = Namespace('http://identifiers.org/hgnc.symbol/')
    up_ns = Namespace('http://identifiers.org/uniprot/')
    g = Graph()

    isa = rn.term('isa')
    partof = rn.term('partof')


    family_names = set([])
    with open(relations_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_MINIMAL,
                                lineterminator='\r\n')
        for line in csv_reader:
            ns1, id1, rel, ns2, id2 = line
            term1 = make_term(ns1, id1)
            term2 = make_term(ns2, id2)
            if rel in ('isa', 'partof'):
                rel_term = rn.term(rel)
            else:
                raise ValueError("Invalid relation %s" % rel)
            g.add((term1, rel_term, term2))

    with open('../resources/entity_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))

