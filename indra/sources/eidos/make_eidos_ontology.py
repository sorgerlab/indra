import yaml
import requests
from os.path import join, dirname, abspath
from rdflib import Graph, Namespace, Literal

eidos_ont_url = 'https://raw.githubusercontent.com/clulab/eidos/master/' + \
                'src/main/resources/org/clulab/wm/eidos/toy_ontology.yml'

eidos_ns = Namespace('https://github.com/clulab/eidos/wiki/JSON-LD/Grounding#')
indra_ns = 'http://sorger.med.harvard.edu/indra/'
indra_rel_ns = Namespace(indra_ns + 'relations/')
isa = indra_rel_ns.term('isa')

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


def get_term(node, prefix):
    node = node.replace(' ', '_')
    path = prefix + '/' + node if prefix else node
    return eidos_ns.term(path)


def build_relations(G, node, tree, prefix):
    this_term = get_term(node, prefix)
    this_prefix = prefix + '/' + node if prefix else node
    for entry in tree:
        if isinstance(entry, str):
            child = entry
        elif isinstance(entry, dict):
            child = list(entry.keys())[0]
            build_relations(G, child, entry[child], this_prefix)
        child_term = get_term(child, this_prefix)
        rel = (child_term, isa, this_term)
        G.add(rel)


if __name__ == '__main__':
    yml = requests.get(eidos_ont_url).content
    root = yaml.load(yml)
    G = Graph()
    for top_entry in root:
        node = list(top_entry.keys())[0]
        build_relations(G, node, top_entry[node], None)
    rdf_path = join(dirname(abspath(__file__)), 'eidos_ontology.rdf')
    save_hierarchy(G, rdf_path)
