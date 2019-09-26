"""This script loads the ontologies for Eidos and Hume and generates RDFs.

The script can handle any ontology which uses the same format (yaml ontology
following the namespace defined at `eidos_ns`).
"""
import yaml
import argparse
import requests
from os.path import join, dirname, abspath
from rdflib import Graph, Namespace, Literal

from indra.sources import eidos, hume


eidos_ns = Namespace('https://github.com/clulab/eidos/wiki/JSON-LD/Grounding#')
indra_ns = 'http://sorger.med.harvard.edu/indra/'
indra_rel_ns = Namespace(indra_ns + 'relations/')
isa = indra_rel_ns.term('isa')

wm_ont_url = ('https://raw.githubusercontent.com/WorldModelers/'
              'Ontologies/master/wm_metadata.yml')
eidos_ont_url = ('https://raw.githubusercontent.com/clulab/eidos/master/'
                 'src/main/resources/org/clulab/wm/eidos/english/'
                 'ontologies/un_ontology.yml')
hume_ont_url = ('https://raw.githubusercontent.com/BBN-E/Hume/master/'
                'resource/ontologies/open/hume_ontology.yaml')


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


def build_relations(G, node, tree, prefix, add_leaves=True):
    this_term = get_term(node, prefix)
    node = node.replace(' ', '_')
    if prefix is not None:
        prefix = prefix.replace(' ', '_')
    this_prefix = prefix + '/' + node if prefix else node
    for entry in tree:
        if isinstance(entry, str):
            if add_leaves:
                child = entry
            else:
                continue
        elif isinstance(entry, dict):
            if 'OntologyNode' not in entry.keys():
                for child in entry.keys():
                    if child[0] != '_' and child != 'examples' \
                       and isinstance(entry[child], (list, dict)):
                        build_relations(
                            G, child, entry[child], this_prefix, add_leaves)
            else:
                child = entry['name']

        if child[0] != '_' and child != 'examples':
            child_term = get_term(child, this_prefix)
            rel = (child_term, isa, this_term)
            G.add(rel)


def update_ontology(ont_url, rdf_path):
    """Load an ontology formatted like Eidos' from github."""
    yaml_root = load_yaml_from_url(ont_url)
    if ont_url == hume_ont_url:
        G = rdf_graph_from_yaml(yaml_root, add_leaves=False)
    else:
        G = rdf_graph_from_yaml(yaml_root, add_leaves=True)
    save_hierarchy(G, rdf_path)


def rdf_graph_from_yaml(yaml_root, add_leaves):
    """Convert the YAML object into an RDF Graph object."""
    G = Graph()
    for top_entry in yaml_root:
        assert len(top_entry) == 1
        node = list(top_entry.keys())[0]
        build_relations(G, node, top_entry[node], None, add_leaves)
    return G


def load_yaml_from_url(ont_url):
    """Return a YAML object loaded from a YAML file URL."""
    res = requests.get(ont_url)
    if res.status_code != 200:
        raise Exception('Could not load ontology from %s' % ont_url)
    root = yaml.load(res.content)
    return root


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process hume or eidos '
        'ontologies. With arguments, only the ontology at <url> is processed '
        'and saved to <fname>. Otherwise the script will process the default '
        'eidos and hume ontologies.')
    parser.add_argument('--url', help='Specify a url to download an ontology '
                                      'from.')
    parser.add_argument('--fname', help='Filename to save new ontology '
                                        'mapping from <url>')
    args = parser.parse_args()

    # Update with from given url to given fname
    if args.fname and args.url:
        update_ontology(args.url, args.fname)

    # Specifying only a URL or only a filename is ambiguous, don't execute
    elif bool(args.fname) ^ bool(args.url):
        print('Must specify both --fname and --url or run without arguments.')

    # Default script execution, backwards compatible
    else:
        # Eidos
        eidos_rdf_path = join(dirname(abspath(eidos.__file__)),
                              'eidos_ontology.rdf')
        update_ontology(eidos_ont_url, eidos_rdf_path)

        # Hume
        hume_rdf_path = join(dirname(abspath(hume.__file__)),
                             'hume_ontology.rdf')
        update_ontology(hume_ont_url, hume_rdf_path)
