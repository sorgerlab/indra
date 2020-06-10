import requests
from .ontology_graph import IndraOntology


class VirtualOntology(IndraOntology):
    def __init__(self, url, ontology='bio'):
        super().__init__()
        self.url = url
        self.ontology = ontology

    def initialize(self):
        self._initialized = True

    def _rel(self, ns, id, rel_types, direction):
        url = self.url + '/%s_rel' % direction
        res = requests.get(url,
                           json={'ns': ns,
                                 'id': id,
                                 'rel_types': rel_types,
                                 'ontology':  self.ontology})
        return res.json()

    def child_rel(self, ns, id, rel_types):
        return self._rel(ns, id, rel_types, 'child')

    def parent_rel(self, ns, id, rel_types):
        return self._rel(ns, id, rel_types, 'parent')

    def get_node_property(self, ns, id, property):
        url = self.url + '/get_node_property'
        res = requests.get(url,
                           json={'ns': ns,
                                 'id': id,
                                 'property': property,
                                 'ontology':  self.ontology})
        return res.json()

