from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import abspath, dirname, join, exists
import rdflib
from rdflib import Namespace, Literal
import pickle

prefixes = """
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX go: <http://purl.obolibrary.org/obo/go#>
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    """

def save_hierarchy(g, path):
    with open(path, 'wb') as out_file:
        g_bytes = g.serialize(format='nt')
        # Replace extra new lines in string and get rid of empty line at end
        g_bytes = g_bytes.replace(b'\n\n', b'\n').strip()
        # Split into rows and sort
        rows = g_bytes.split(b'\n')
        rows.sort()
        g_bytes = b'\n'.join(rows)
        out_file.write(g_bytes)

def get_cellular_components(g):
    # Query for direct part_of relationships
    query = prefixes + """
        SELECT ?id ?label ?supid ?suplabel
        WHERE {
            ?class oboInOwl:hasOBONamespace "cellular_component"^^xsd:string .
            ?class oboInOwl:id ?id .
            ?class rdfs:label ?label .
            ?class rdfs:subClassOf ?sup .
            ?sup oboInOwl:hasOBONamespace "cellular_component"^^xsd:string .
            ?sup oboInOwl:id ?supid .
            ?sup rdfs:label ?suplabel
            }
        """
    res1 = g.query(query)
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
            ?sup oboInOwl:hasOBONamespace "cellular_component"^^xsd:string .
            ?sup oboInOwl:id ?supid .
            ?sup rdfs:label ?suplabel
            }
        """
    res2 = g.query(query)
    res = list(res1) + list(res2)
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
    ln = Namespace(indra_ns + 'locations/')
    rn = Namespace(indra_ns + 'relations/')
    part_of = rn.term('partof')
    has_name = rn.term('hasName')
    for comp_id, comp_name in component_map.items():
        g.add((ln.term(comp_id), has_name, Literal(comp_name)))
        sups = component_part_map.get(comp_id)
        if sups is not None:
            for sup_id in sups:
                g.add((ln.term(comp_id), part_of, ln.term(sup_id)))
    return g


def main():
    # This file can be donwloaded from:
    # http://geneontology.org/ontology/go.owl
    go_owl_file = join(dirname(abspath(__file__)), '../../data/go.owl')
    pkl_file = join(dirname(abspath(__file__)), '../../data/go.pkl')
    rdf_file = join(dirname(abspath(__file__)),
                    '../resources/cellular_component_hierarchy.rdf')
    g = rdflib.Graph()
    if not exists(pkl_file):
        print('Parsing %s' % go_owl_file)
        g.parse(abspath(go_owl_file))
        with open(pkl_file, 'wb') as fh:
            pickle.dump(g, fh, protocol=2)
    else:
        print('Loading %s' % pkl_file)
        with open(pkl_file, 'rb') as fh:
            g = pickle.load(fh)
    print('Getting cellular components')
    component_map, component_part_map = get_cellular_components(g)
    gg = make_component_hierarchy(component_map, component_part_map)
    save_hierarchy(gg, rdf_file)

if __name__ == '__main__':
    main()
