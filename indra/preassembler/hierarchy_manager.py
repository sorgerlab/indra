from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import rdflib
import logging
try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache

from indra.preassembler.make_entity_hierarchy import ns_map

logger = logging.getLogger('hierarchy_manager')

class HierarchyManager(object):
    """Store hierarchical relationships between different types of entities.

    Used to store, e.g., entity hierarchies (proteins and protein families)
    and modification hierarchies (serine phosphorylation vs. phosphorylation).

    Parameters
    ----------
    rdf_file : string
        Path to the RDF file containing the hierarchy.
    build_closure : Optional[bool]
        If True, the transitive closure of the hierarchy is generated
        up from to speed up processing. Default: True
    uri_as_name: Optional[bool]
        If True, entries are accessed directly by their URIs. If False
        entries are accessed by finding their name through the
        hasName relationship. Default: True

    Attributes
    ----------
    graph : instance of `rdflib.Graph`
        The RDF graph containing the hierarchy.
    """
    prefixes = """
        PREFIX rn: <http://sorger.med.harvard.edu/indra/relations/>
        """

    def __init__(self, rdf_file, build_closure=True, uri_as_name=True):
        """Initialize with the path to an RDF file"""
        self.build_closure = build_closure
        self.uri_as_name = uri_as_name
        self.graph = rdflib.Graph()
        self.graph.parse(rdf_file, format='nt')
        self.relations_prefix = \
            'http://sorger.med.harvard.edu/indra/relations/'
        self.initialize()


    def initialize(self):
        self.isa_closure = {}
        self.partof_closure = {}
        self.components = {}
        if self.build_closure:
            self.build_transitive_closures()
        # Build reverse lookup dict from the entity hierarchy
        self._children = {}
        all_children = set(self.isa_closure.keys()).union(
                            self.partof_closure.keys())
        for child in all_children:
            parents = self.get_parents(child)
            for parent in parents:
                children_list = self._children.get(parent, [])
                children_list.append(child)
                self._children[parent] = children_list

    def build_transitive_closures(self):
        """Build the transitive closures of the hierarchy.

        This method constructs dictionaries which contain terms in the
        hierarchy as keys and either all the "isa+" or "partof+" related terms
        as values.
        """
        self.component_counter = 0
        for rel, tc_dict in (('isa', self.isa_closure),
                             ('partof', self.partof_closure)):
            rel_uri = self.relations_prefix + rel
            rel_ref = rdflib.term.URIRef(rel_uri)
            for x in self.graph.all_nodes():
                rel_closure = self.graph.transitive_objects(x, rel_ref)
                xs = x.toPython()
                for y in rel_closure:
                    ys = y.toPython()
                    if xs == ys:
                        continue
                    try:
                        tc_dict[xs].append(ys)
                    except KeyError:
                        tc_dict[xs] = [ys]
                    self._add_component(xs, ys)

    def _add_component(self, xs, ys):
        xcomp = self.components.get(xs)
        ycomp = self.components.get(ys)
        if xcomp is None:
            if ycomp is None:
                # Neither x nor y are in a component so we start a
                # new component and assign x and y to the same
                # component
                self.components[xs] = self.component_counter
                self.components[ys] = self.component_counter
                self.component_counter += 1
            else:
                # Because y is already part of an existing component
                # we assign its component to x
                self.components[xs] = ycomp
        else:
            if ycomp is None:
                # Because x is already part of an existing component
                # we assign its component to y
                self.components[ys] = xcomp
            else:
                # This is a special case in which both x and y are
                # parts of components
                # If they are in the same component then there's
                # nothing further to do
                if xcomp != ycomp:
                    remove_component = max(xcomp, ycomp)
                    joint_component = min(xcomp, ycomp)
                    for k, v in self.components.items():
                        if v == remove_component:
                            self.components[k] = joint_component


    @lru_cache(maxsize=100000)
    def find_entity(self, x):
        """
        Get the entity that has the specified name (or synonym).

        Parameters
        ----------
        x : string
            Name or synonym for the target entity.
        """

        qstr = self.prefixes + """
            SELECT ?x WHERE {{
                ?x rn:hasName "{0}" .
            }}
            """.format(x)
        res = self.graph.query(qstr)
        if list(res):
            en = list(res)[0][0].toPython()
            return en
        else:
            return None


    def isa(self, ns1, id1, ns2, id2):
        """Indicate whether one entity has an "isa" relationship to another.

        Parameters
        ----------
        ns1 : string
            Namespace code for an entity.
        id1 : string
            URI for an entity.
        ns2 : string
            Namespace code for an entity.
        id2 : string
            URI for an entity.

        Returns
        -------
        bool
            True if t1 has an "isa" relationship with t2, either directly or
            through a series of intermediates; False otherwise.
        """
        # if id2 is None, or both are None, then it's by definition isa:
        if id2 is None or (id2 is None and id1 is None):
            return True
        # If only id1 is None, then it cannot be isa
        elif id1 is None:
            return False

        if self.isa_closure:
            term1 = self.get_uri(ns1, id1)
            term2 = self.get_uri(ns2, id2)
            ec = self.isa_closure.get(term1)
            if ec is not None and term2 in ec:
                return True
            else:
                return False
        else:
            if not self.uri_as_name:
                e1 = self.find_entity(id1)
                e2 = self.find_entity(id2)
                if e1 is None or e2 is None:
                    return False
                t1 = rdflib.term.URIRef(e1)
                t2 = rdflib.term.URIRef(e2)
            else:
                u1 = self.get_uri(ns1, id1)
                u2 = self.get_uri(ns2, id2)
                t1 = rdflib.term.URIRef(u1)
                t2 = rdflib.term.URIRef(u2)

            rel_uri = self.relations_prefix + 'isa'
            rel_ref = rdflib.term.URIRef(rel_uri)
            to = self.graph.transitive_objects(t1, rel_ref)
            if t2 in to:
                return True
            else:
                return False

    def partof(self, ns1, id1, ns2, id2):
        """Indicate whether one entity is physically part of another.

        Parameters
        ----------
        ns1 : string
            Namespace code for an entity.
        id1 : string
            URI for an entity.
        ns2 : string
            Namespace code for an entity.
        id2 : string
            URI for an entity.

        Returns
        -------
        bool
            True if t1 has a "partof" relationship with t2, either directly or
            through a series of intermediates; False otherwise.
        """
        # if id2 is None, or both are None, then it's by definition isa:
        if id2 is None or (id2 is None and id1 is None):
            return True
        # If only id1 is None, then it cannot be isa
        elif id1 is None:
            return False

        if self.partof_closure:
            term1 = self.get_uri(ns1, id1)
            term2 = self.get_uri(ns2, id2)
            ec = self.partof_closure.get(term1)
            if ec is not None and term2 in ec:
                return True
            else:
                return False
        else:
            if not self.uri_as_name:
                e1 = self.find_entity(id1)
                e2 = self.find_entity(id2)
                if e1 is None or e2 is None:
                    return False
                t1 = rdflib.term.URIRef(e1)
                t2 = rdflib.term.URIRef(e2)
            else:
                u1 = self.get_uri(ns1, id1)
                u2 = self.get_uri(ns2, id2)
                t1 = rdflib.term.URIRef(u1)
                t2 = rdflib.term.URIRef(u2)

            rel_uri = self.relations_prefix + 'partof'
            rel_ref = rdflib.term.URIRef(rel_uri)
            to = self.graph.transitive_objects(t1, rel_ref)
            if t2 in to:
                return True
            else:
                return False

    def get_parents(self, uri, type='all'):
        """Return parents of a given entry.

        Parameters
        ----------
        uri : str
            The URI of the entry whose parents are to be returned. See the
            get_uri method to construct this URI from a name space and id.
        type : str
            'all': return all parents irrespective of level;
            'immediate': return only the immediate parents;
            'top': return only the highest level parents
        """
        immediate_parents = set(self.isa_closure.get(uri, [])).union(
                                set(self.partof_closure.get(uri, [])))
        if type == 'immediate':
            return immediate_parents
        all_parents = set()
        for parent in immediate_parents:
            grandparents = self.get_parents(parent, type='all')
            all_parents = all_parents.union(grandparents)
        all_parents = all_parents.union(immediate_parents)
        if type == 'all':
            return all_parents
        else:
            top_parents = set()
            for parent in all_parents:
                if not self.get_parents(parent, type='immediate'):
                    top_parents.add(parent)
            return top_parents
        return

    def get_children(self, uri):
        """Return all (not just immediate) children of a given entry.

        Parameters
        ----------
        uri : str
            The URI of the entry whose children are to be returned. See the
            get_uri method to construct this URI from a name space and id.
        """
        children = self._children.get(uri, [])
        return children

    @lru_cache(maxsize=100000)
    def query_rdf(self, id1, rel, id2):
        term1 = self.find_entity(id1)
        term2 = self.find_entity(id2)
        qstr = self.prefixes + """ 
            SELECT (COUNT(*) as ?s) WHERE {{
                <{}> {} <{}> .
                }}
            """.format(term1, rel, term2)
        res = self.graph.query(qstr)
        count = [r[0] for r in res][0]
        if count.toPython() == 1:
            return True
        else:
            return False

    @staticmethod
    def get_uri(ns, id):
        if ns == 'HGNC':
            return 'http://identifiers.org/hgnc.symbol/' + id
        elif ns == 'UP':
            return 'http://identifiers.org/uniprot/' + id
        elif ns == 'FPLX' or ns == 'INDRA':
            return 'http://sorger.med.harvard.edu/indra/entities/' + id
        else:
            return ns + id

    @staticmethod
    def ns_id_from_uri(uri):
        sep_ix = uri.rfind('/') + 1
        ag_ns = uri[0:sep_ix]
        ag_id = uri[sep_ix:]
        # Handle one special case here for HGNC IDs
        if ag_id.startswith('HGNC:'):
            ag_ns = 'http://identifiers.org/hgnc.symbol/'
            ag_id = hgnc_client.get_hgnc_name(db_id[5:])
        ag_ns_name = ns_map.get(ag_ns)
        if ag_ns_name is None:
            raise UnknownNamespaceException('Unknown namespace %s' % ag_ns)
        return (ag_ns_name, ag_id)

# Load the default entity and modification hierarchies
entity_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/entity_hierarchy.rdf')
mod_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/modification_hierarchy.rdf')
act_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/activity_hierarchy.rdf')
ccomp_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/cellular_component_hierarchy.rdf')

# Default entity hierarchy loaded from the RDF file at
# `resources/entity_hierarchy.rdf`.
entity_hierarchy = HierarchyManager(entity_file_path, build_closure=True,
                                    uri_as_name=True)
# Default modification hierarchy loaded from the RDF file at
# `resources/modification_hierarchy.rdf`.
modification_hierarchy = HierarchyManager(mod_file_path, build_closure=True,
                                          uri_as_name=True)
# Default activity hierarchy loaded from the RDF file at
# `resources/activity_hierarchy.rdf`.
activity_hierarchy = HierarchyManager(act_file_path, build_closure=True,
                                      uri_as_name=True)
# Default cellular_component hierarchy loaded from the RDF file at
# `resources/cellular_component_hierarchy.rdf`.
ccomp_hierarchy = HierarchyManager(ccomp_file_path, build_closure=False,
                                   uri_as_name=False)

hierarchies = {'entity': entity_hierarchy,
               'modification': modification_hierarchy,
               'activity': activity_hierarchy,
               'cellular_component': ccomp_hierarchy}
