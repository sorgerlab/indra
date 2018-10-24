import os
import rdflib

fname = '../../data/go.owl'

g = rdflib.Graph()
print("Parsing graph")
g.parse(os.path.abspath(fname))

prefixes = """
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX go: <http://purl.obolibrary.org/obo/go#>
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    """

query = prefixes + """
    SELECT ?id ?label
    WHERE {
        ?class oboInOwl:id ?id .
        ?class rdfs:label ?label 
    }
"""
print("Running query")
res = g.query(query)
mappings = []
for id_lit, label_lit in res:
    mappings.append((id_lit.value, label_lit.value))

