"""This script loads the ontologies for Eidos and Hume and generates RDFs.

The script can handle any ontology which uses the same format (yaml ontology
following the namespace defined at `eidos_ns`).
"""
import os
import yaml
import logging
import argparse
import requests
from rdflib import Graph, Namespace, Literal


logger = logging.getLogger('indra.preassembler.make_wm_ontologies')


eidos_ns = Namespace('https://github.com/clulab/eidos/wiki/JSON-LD/Grounding#')
indra_ns = 'http://sorger.med.harvard.edu/indra/'
indra_rel_ns = Namespace(indra_ns + 'relations/')
isa = indra_rel_ns.term('isa')
isequal = indra_rel_ns.term('is_equal')
isopp = indra_rel_ns.term('is_opposite')

wm_ont_url = ('https://raw.githubusercontent.com/WorldModelers/'
              'Ontologies/master/wm_metadata.yml')
eidos_ont_url = ('https://raw.githubusercontent.com/clulab/eidos/master/'
                 'src/main/resources/org/clulab/wm/eidos/english/'
                 'ontologies/un_ontology.yml')
hume_ont_url = ('https://raw.githubusercontent.com/BBN-E/Hume/master/'
                'resource/ontologies/open/hume_ontology.yaml')


class HierarchyConverter(object):
    def __init__(self, url, path, add_leaves=True):
        self.yml = load_yaml_from_url(url)
        self.path = path
        self.add_leaves = add_leaves
        self.G = None

    def convert_ontology(self):
        self.G = Graph()
        for top_entry in self.yml:
            node = list(top_entry.keys())[0]
            self.build_relations(node, top_entry[node], None)

    def build_relations(self, node, tree, prefix):
        this_term = get_term(node, prefix)
        node = node.replace(' ', '_')
        if prefix is not None:
            prefix = prefix.replace(' ', '_')
        this_prefix = prefix + '/' + node if prefix else node
        for entry in tree:
            if isinstance(entry, str):
                if self.add_leaves:
                    child = entry
                else:
                    continue
            elif isinstance(entry, dict):
                if 'OntologyNode' not in entry.keys():
                    for child in entry.keys():
                        if child[0] != '_' and child != 'examples' \
                                and isinstance(entry[child], (list, dict)):
                            self.build_relations(child, entry[child],
                                                 this_prefix)
                else:
                    child = entry['name']

            if child[0] != '_' and child != 'examples':
                child_term = get_term(child, this_prefix)
                rel = (child_term, isa, this_term)
                self.G.add(rel)
                opp = entry.get('opposite')
                if opp:
                    parts = opp.split('/')
                    opp_term = get_term(parts[-1], '/'.join(parts[:-1]))
                    rel = (opp_term, isopp, child_term)
                    self.G.add(rel)

    def save_hierarchy(self):
        g_bytes = self.G.serialize(format='nt')
        # Replace extra new lines in string and get rid of empty
        # line at end
        g_bytes = g_bytes.replace(b'\n\n', b'\n').strip()
        # Split into rows and sort
        rows = g_bytes.split(b'\n')
        rows.sort()
        g_bytes = b'\n'.join(rows)
        with open(self.path, 'wb') as out_file:
            out_file.write(g_bytes)


def get_term(node, prefix):
    node = node.replace(' ', '_')
    path = prefix + '/' + node if prefix else node
    return eidos_ns.term(path)


def load_yaml_from_url(ont_url):
    """Return a YAML object loaded from a YAML file URL."""
    res = requests.get(ont_url)
    if res.status_code != 200:
        raise Exception('Could not load ontology from %s' % ont_url)
    root = yaml.load(res.content, Loader=yaml.FullLoader)
    return root


if __name__ == '__main__':
    wm_rdf_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               os.pardir, 'resources', 'wm_ontology.rdf')
    parser = argparse.ArgumentParser(
        description='Process YAML ontology files to create INDRA-compatible '
                    'RDF files as input to the Preassembler.')
    parser.add_argument('--url', help='Specify the url to download an '
                                      'ontology from.',
                        default=wm_ont_url)
    parser.add_argument('--fname', help='Name of the file to save new RDF '
                                        'ontology into.',
                        default=wm_rdf_path)
    args = parser.parse_args()

    hc = HierarchyConverter(args.url, args.fname, args.url == hume_ont_url)
    hc.convert_ontology()
    hc.save_hierarchy()
