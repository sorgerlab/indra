import os
import rdflib
from rdflib import Namespace, Literal
import cPickle as pickle

prefixes = """
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX go: <http://purl.obolibrary.org/obo/go#>
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    """

class CellularComponent(object):
    def __init__(self, go_id, name):
        self.go_id = go_id
        self.name = name

def get_cellular_components(g):
    # Query for direct part_of relationships
    query = prefixes + """
        SELECT ?id ?label ?supid ?suplabel
        WHERE {
            ?class oboInOwl:hasOBONamespace "cellular_component"^^xsd:string .
            ?class oboInOwl:id ?id .
            ?class rdfs:label ?label .
            ?class rdfs:subClassOf ?restr .
            ?restr owl:onProperty ?prop .
            ?prop oboInOwl:id "part_of"^^xsd:string .
            ?restr owl:someValuesFrom ?sup .
            ?sup oboInOwl:id ?supid .
            ?sup rdfs:label ?suplabel
            }
        """
    res = g.query(query)
    component_map = {}
    component_part_map = {}
    for r in res:
        comp_id, comp_name, sup_id, sup_name = [rr.toPython() for rr in r]
        component_map[comp_id] = comp_name
        component_map[sup_id] = sup_name
        try:
            component_part_map[comp_id].append(sup_id)
        except KeyError:
            component_part_map[comp_id] = [sup_id]
    # Query for isa + part_of relationships
    query = prefixes + """
        SELECT ?id ?label ?supid ?suplabel
        WHERE {
            ?class oboInOwl:hasOBONamespace "cellular_component"^^xsd:string .
            ?class oboInOwl:id ?id .
            ?class rdfs:label ?label .
            ?class rdfs:subClassOf+ ?supclass .
            ?supclass oboInOwl:hasOBONamespace "cellular_component"^^xsd:string .
            ?supclass rdfs:subClassOf ?restr .
            ?restr owl:onProperty ?prop .
            ?prop oboInOwl:id "part_of"^^xsd:string .
            ?restr owl:someValuesFrom ?sup .
            ?sup oboInOwl:id ?supid .
            ?sup rdfs:label ?suplabel
            }
        """
    res = g.query(query)
    component_map = {}
    component_part_map = {}
    for r in res:
        comp_id, comp_name, sup_id, sup_name = [rr.toPython() for rr in r]
        component_map[comp_id] = comp_name
        component_map[sup_id] = sup_name
        try:
            component_part_map[comp_id].append(sup_id)
        except KeyError:
            component_part_map[comp_id] = [sup_id]
    return component_map, component_part_map

def make_component_hierarchy(component_map, component_part_map):
    g = rdflib.Graph()
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    en = Namespace(indra_ns + 'entities/')
    rn = Namespace(indra_ns + 'relations/')
    part_of = rn.term('partof')
    has_name = rn.term('hasName')
    for comp_id, comp_name in component_map.iteritems():
        g.add((en.term(comp_id), has_name, Literal(comp_name)))
        sups = component_part_map.get(comp_id)
        if sups is not None:
            for sup_id in sups:
                g.add((en.term(comp_id), part_of, en.term(sup_id)))
    return g

if __name__ == '__main__':
    # This file can be donwloaded from:
    # http://purl.obolibrary.org/obo/go.owl
    fname = '/home/beni/data/go.owl'
    g = rdflib.Graph()
    pkl_name = '/home/beni/data/go.pkl'
    if not os.path.exists(pkl_name):
        print 'Parsing %s' % fname
        g.parse(fname)
        with open(pkl_name, 'wb') as fh:
            pickle.dump(g, fh)
    else:
        print 'Loading %s' % pkl_name
        g = pickle.load(open(pkl_name, 'rb'))
    print 'Getting cellular components'
    component_map, component_part_map = get_cellular_components(g)
    gg = make_component_hierarchy(component_map, component_part_map)
    with open('cellular_component_hierarchy.rdf', 'wt') as out_file:
        out_file.write(gg.serialize(format='xml'))
