import logging
import networkx
from collections import deque
from .shortest_path import bidirectional_shortest_path

logger = logging.getLogger(__name__)


class IndraOntology(networkx.MultiDiGraph):
    def _check_path(self, ns1, id1, ns2, id2, edge_types):
        try:
            path = bidirectional_shortest_path(self,
                                               label(ns1, id1),
                                               label(ns2, id2),
                                               edge_types=edge_types)
        except networkx.NetworkXError:
            return False
        return True

    @staticmethod
    def get_ns_id(node):
        return node.split(':', maxsplit=1)

    @staticmethod
    def get_ns(node):
        return IndraOntology.get_ns_id(node)[0]

    @staticmethod
    def get_id(node):
        return IndraOntology.get_ns_id(node)[1]

    def isrel(self, ns1, id1, ns2, id2, rels):
        return self._check_path(ns1, id1, ns2, id2, rels)

    def isa(self, ns1, id1, ns2, id2):
        return self.isrel(ns1, id1, ns2, id2, rels={'isa'})

    def partof(self, ns1, id1, ns2, id2):
        return self.isrel(ns1, id1, ns2, id2, rels={'partof'})

    def isa_or_partof(self, ns1, id1, ns2, id2):
        return self.isrel(ns1, id1, ns2, id2, rels={'isa', 'partof'})

    def maps_to(self, ns1, id1, ns2, id2):
        return self._check_path(ns1, id1, ns2, id2, {'xref'})

    def map_to(self, ns1, id1, ns2):
        source = label(ns1, id1)
        targets = [target for target in self.xrefs(source)
                   if self.get_ns(target) == ns2]
        if len(targets) == 1:
            return targets[0]
        return None

    def xrefs(self, source):
        for _, target, data in self.edges(source, data=True):
            if data['type'] == 'xref':
                yield target

    def get_mappings(self, ns, id):
        source = label(ns, id)
        visited = {source}
        queue = deque([(source, self.xrefs(source))])
        targets = []
        while queue:
            parent, children = queue[0]
            try:
                child = next(children)
                if child not in visited:
                    targets.append(child)
                    visited.add(child)
                    queue.append((child, self.xrefs(child)))
            except StopIteration:
                queue.popleft()
        return [self.get_ns_id(t) for t in targets]

    def get_name(self, ns, id):
        return self.get_node_property(ns, id, property='name')

    def get_polarity(self, ns, id):
        return self.get_node_property(ns, id, property='polarity')

    def get_node_property(self, ns, id, property):
        try:
            return self.nodes[label(ns, id)][property]
        except KeyError:
            return None


def label(ns, id):
    return '%s:%s' % (ns, id)


