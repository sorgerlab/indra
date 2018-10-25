import os
import csv
from os.path import abspath, dirname, join
import rdflib
from indra.util import read_unicode_csv, write_unicode_csv


go_file = join(dirname(abspath(__file__)), '..', 'resources',
                 'go_id_label_mappings.tsv')


go_mappings = {}
for go_id, go_label in read_unicode_csv(go_file, delimiter='\t'):
    go_mappings[go_id] = go_label


def get_mappings_from_owl():
    """Compile all ID->label mappings from the GO OWL file."""
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
    # Write to file
    write_unicode_csv(go_file, mappings)


def get_go_label(go_id):
    return go_mappings.get(go_id)

