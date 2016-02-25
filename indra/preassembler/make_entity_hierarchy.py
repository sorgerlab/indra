import sys
from rdflib import Graph, Namespace, Literal
import csv
import urllib2

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
    has_long_name = rn.term('hasLongName')
    has_synonym = rn.term('hasSynonym')
    isa = rn.term('isa')

    # Read BEL family names
    res = urllib2.urlopen('http://resource.belframework.org/belframework/'+\
        'latest-release/namespace/selventa-protein-families.belns')
    belns_text = res.read()
    start = belns_text.find('[Values]')
    lines = belns_text[start:].split('\n')
    bel_family_names = []
    for l in lines:
        if l.endswith(' Family|P'):
            family_name = l[:-2].replace(' ', '_').replace('/', '_')
            bel_family_names.append(family_name)

    family_names = set([])
    with open(proteins_file) as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            print line
            symbol, name, synonyms_str, family = line
            g.add((en.term(symbol), has_long_name, Literal(name.strip())))
            g.add((en.term(symbol), has_name, Literal(symbol)))
            synonyms = synonyms_str.split(', ')
            for s in synonyms:
                if s != '':
                    g.add((en.term(symbol), has_synonym, Literal(s)))
            family_name = family.strip()
            if family_name != '':
                print family_name
                g.add((en.term(symbol), isa, en.term(family_name)))
                g.add((en.term(family), has_name, Literal(family_name)))
                family_names.add(family_name)
    for fn in family_names:
        bel_synonyms = [bn for bn in bel_family_names 
                        if bn.find(fn) != -1]
        for b in bel_synonyms:
            g.add((en.term(fn), has_synonym, Literal(b.upper())))

    # Further BEL family names added manually
    g.add((en.term('ERK'), has_synonym, Literal('MAPK_ERK1_2_FAMILY')))

    with open('entity_hierarchy.rdf', 'wt') as out_file:
        out_file.write(g.serialize(format='xml'))

