"""This script loads the ontologies for Eidos and Hume and generates RDFs.

The script can handle any ontology which uses the same format (yaml ontology
following the namespace defined at `eidos_ns`).
"""
import logging
import requests
from collections import defaultdict
from indra.pipeline import register_pipeline
from ..ontology_graph import IndraOntology, with_initialize


logger = logging.getLogger(__name__)

wm_ont_url = ('https://raw.githubusercontent.com/WorldModelers/'
              'Ontologies/master/wm_flat_metadata.yml')


def get_term(node, prefix):
    node = node.replace(' ', '_')
    path = prefix + '/' + node if prefix else node
    return WorldOntology.label('WM', path)


def load_yaml_from_url(ont_url):
    """Return a YAML object loaded from a YAML file URL."""
    import yaml
    res = requests.get(ont_url)
    res.raise_for_status()
    root = yaml.load(res.content, Loader=yaml.FullLoader)
    return root


class WorldOntology(IndraOntology):
    """Represents the ontology used for World Modelers applications.

    Parameters
    ----------
    url : str
        The URL pointing to a World Modelers ontology YAML.

    Attributes
    ----------
    url : str
        The URL pointing to a World Modelers ontology YAML.
    yml : list
        The ontology YAML as loaded by the yaml package from the
        URL.
    """
    name = 'world'
    version = '1.0'

    def __init__(self, url):
        super().__init__()
        self.yml = None
        self.url = url

    def initialize(self):
        """Load the World Modelers ontology from the web and build the
        graph."""
        logger.info('Initializing world ontology from %s' % self.url)
        self.add_wm_ontology(self.url)
        self._initialized = True
        self._build_transitive_closure()
        logger.info('Ontology has %d nodes' % len(self))

    def add_wm_ontology(self, url):
        self.yml = load_yaml_from_url(url)
        self._load_yml(self.yml)

    @with_initialize
    def dump_yml_str(self):
        """Return a string-serialized form of the loaded YAML

        Returns
        -------
        str
            The YAML string of the ontology.
        """
        import yaml
        return yaml.dump(self.yml)

    def _load_yml(self, yml):
        self.clear()
        for top_entry in yml:
            node = list(top_entry.keys())[0]
            self.build_relations(node, top_entry[node], None)

    def build_relations(self, node, tree, prefix):
        nodes = defaultdict(dict)
        edges = []
        this_term = get_term(node, prefix)
        node = node.replace(' ', '_')
        if prefix is not None:
            prefix = prefix.replace(' ', '_')
        this_prefix = prefix + '/' + node if prefix else node
        for entry in tree:
            if isinstance(entry, str):
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
                edges.append((child_term, this_term, {'type': 'isa'}))
                opp = entry.get('opposite')
                if opp:
                    parts = opp.split('/')
                    opp_term = get_term(parts[-1], '/'.join(parts[:-1]))
                    edges.append((opp_term, child_term, {'type': 'is_opposite'}))
                    edges.append((child_term, opp_term, {'type': 'is_opposite'}))
                pol = entry.get('polarity')
                if pol is not None:
                    nodes[child_term]['polarity'] = pol
        self.add_nodes_from([(k, v) for k, v in dict(nodes).items()])
        self.add_edges_from(edges)

    def add_entry(self, entry, examples=None):
        """Add a new ontology entry with examples.

        This works by adding the entry to the yml attribute first
        and then reloading the entire yaml to build a new graph.

        Parameters
        ----------
        entry : str
            The new entry.
        examples : list of str
            Examples for the new entry.
        """
        examples = examples if examples else []
        parts = entry.split('/')
        root = self.yml
        for idx, part in enumerate(parts):
            new_root = None
            for element in self.yml:
                # If this is an OntologyNode
                if 'OntologyNode' in element:
                    if element['name'] == part:
                        new_root = element
                        break
                else:
                    assert len(element) == 1
                    key = list(element.keys())[0]
                    if key == part:
                        new_root = element[key]
                        break
            if new_root is None:
                if idx == len(parts) - 1:
                    root.append({'OntologyNode': None, 'name': part,
                                 'examples': examples})
                    break
                else:
                    root.append({part: []})
                    new_root = root[-1][part]
            root = new_root
        self._load_yml(self.yml)


@register_pipeline
def load_world_ontology(url=wm_ont_url):
    return WorldOntology(url)


world_ontology = load_world_ontology()
