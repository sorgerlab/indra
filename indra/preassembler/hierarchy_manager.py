import os
import rdflib
import logging
from functools import lru_cache


from indra.preassembler.make_entity_hierarchy import ns_map

logger = logging.getLogger(__name__)


class HierarchyManager(object):
    """Store hierarchical relationships between different types of entities.

    Used to store, e.g., entity hierarchies (proteins and protein families)
    and modification hierarchies (serine phosphorylation vs. phosphorylation).

    Parameters
    ----------
    rdf_file : string
        Path to the RDF file containing the hierarchy.
    build_closure : bool or list or None
        If True, the transitive closure of the hierarchy is generated
        up front to speed up processing. If a list, the entries in the list
        are namespaces for which a transitive closure should be built.
        Namespaces not listed are ignored and isa/partof lookups will not use
        the transitive closure. If False, no transitive closure is built.
        Default: True
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

    def __init__(self, rdf_file=None, build_closure=None, uri_as_name=True):
        """Initialize with the path to an RDF file"""
        self.build_closure = build_closure
        self.uri_as_name = uri_as_name
        self.relations_prefix = \
            'http://sorger.med.harvard.edu/indra/relations/'
        self.isa_closure = set()
        self.partof_closure = set()
        self.isa_or_partof_closure = set()
        self.components = {}
        self._children = {}
        self.component_counter = 0
        # If an RDF file was given, we build up the internal data structures.
        # Otherwise we defer initialization until later.
        if rdf_file:
            self.load_from_rdf_file(rdf_file)
        else:
            self.graph = None

    def load_from_rdf_file(self, rdf_file):
        """Initialize given an RDF input file representing the hierarchy."

        Parameters
        ----------
        rdf_file : str
            Path to an RDF file.
        """
        self.graph = rdflib.Graph()
        logger.info('Loading hierarchy from %s' % rdf_file)
        self.graph.parse(os.path.abspath(rdf_file), format='nt')
        logger.info('Loaded hierarchy from %s' % rdf_file)
        self.initialize()
        logger.info('Initialized hierarchy')

    def load_from_rdf_string(self, rdf_str):
        """Initialize given an RDF string representing the hierarchy."

        Parameters
        ----------
        rdf_str : str
            An RDF string.
        """
        self.graph = rdflib.Graph()
        self.graph.parse(data=rdf_str, format='nt')
        self.initialize()

    def load_from_rdf_graph(self, rdf_graph):
        """Initialize given an RDF Graph representing the hierarchy."

        Parameters
        ----------
        rdf_graph : rdflib.Graph
            An rdflib Graph representing the hierarchy.
        """
        self.graph = rdf_graph
        self.initialize()

    def initialize(self):
        self.build_transitive_closures()

        # Build reverse lookup dict from the hierarchy
        # First get all URIs that correspond to parents
        all_parents = {parent for child, parent in self.isa_or_partof_closure}
        # We use the inverse relation here
        rel_fun = lambda node, graph: self.isa_or_partof_objects(node,
                                                                 inverse=True)
        # Now for each parent we get the inverse transitive closure to
        # get all its children nodes
        self._children = {}
        for parent in all_parents:
            children = self.graph.transitiveClosure(rel_fun,
                                                    rdflib.term.URIRef(parent))
            children_uris = list(set(c.toPython() for c in children))
            if children_uris:
                self._children[parent] = children_uris

    def extend_with(self, rdf_file):
        """Extend the RDF graph of this HierarchyManager with another RDF file.

        Parameters
        ----------
        rdf_file : str
            An RDF file which is parsed such that the current graph and the
            graph described by the file are merged.
        """
        self.graph.parse(os.path.abspath(rdf_file), format='nt')
        self.initialize()

    def build_transitive_closures(self):
        """Build the transitive closures of the hierarchy.

        This method constructs dictionaries which contain terms in the
        hierarchy as keys and either all the "isa+" or "partof+" related terms
        as values.
        """
        self.component_counter = 0
        for rel, tc_set in ((self.isa_objects, self.isa_closure),
                            (self.partof_objects, self.partof_closure),
                            (self.isa_or_partof_objects,
                             self.isa_or_partof_closure)):
            self.build_transitive_closure(rel, tc_set)

    def build_transitive_closure(self, rel, tc_set):
        """Build a transitive closure for a given relation in a given dict."""
        # If there are no namespaces to build closures for, we
        # return immediately
        if not self.build_closure:
            return
        nodes = self._get_build_transitive_closure_objects()
        # Make a function with the righ argument structure
        rel_fun = lambda node, graph: rel(node)
        for x in nodes:
            rel_closure = self.graph.transitiveClosure(rel_fun, x)
            xs = x.toPython()
            for y in rel_closure:
                ys = y.toPython()
                if xs == ys:
                    continue
                tc_set.add((xs, ys))
                if rel == self.isa_or_partof_objects:
                    self._add_component(xs, ys)

    def _term_in_closure_namespace(self, term):
        """Return True if term is in a namespace with a closure."""
        if self.build_closure is True:
            return True
        elif isinstance(self.build_closure, (list, tuple)):
            return any([term.startswith(ns) for ns in self.build_closure])
        else:
            return False

    def _get_build_transitive_closure_objects(self):
        """Return objects that will be included in the transitive closures."""
        if not self.build_closure:
            nodes = []
        else:
            nodes = list(self.graph.all_nodes())
            if isinstance(self.build_closure, (list, tuple)):
                nodes = [node for node in nodes if
                         self._term_in_closure_namespace(node)]
        return nodes

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

    def isa_objects(self, node, inverse=False):
        # Normally we look for objects of the relation, but if inverted,
        # we look for the subject
        predicate = rdflib.term.URIRef(self.relations_prefix + 'isa')
        partner_list = self.graph.subjects(predicate, node) if inverse else \
            self.graph.objects(node, predicate)
        for o in partner_list:
            yield o

    def partof_objects(self, node, inverse=False):
        # Normally we look for objects of the relation, but if inverted,
        # we look for the subject
        predicate = rdflib.term.URIRef(self.relations_prefix + 'partof')
        partner_list = self.graph.subjects(predicate, node) if inverse else \
            self.graph.objects(node, predicate)
        for o in partner_list:
            yield o

    def isa_or_partof_objects(self, node, inverse=False):
        for o in self.isa_objects(node, inverse):
            yield o
        for o in self.partof_objects(node, inverse):
            yield o

    def directly_or_indirectly_related(self, ns1, id1, ns2, id2, closure_set,
                                       relation_func):
        """Return True if two entities have the speicified relationship.

        This relation is constructed possibly through multiple links connecting
        the two entities directly or indirectly.

        Parameters
        ----------
        ns1 : str
            Namespace code for an entity.
        id1 : str
            URI for an entity.
        ns2 : str
            Namespace code for an entity.
        id2 : str
            URI for an entity.
        closure_set : set
            A set containing tuples of entities that have the
            specified relationship, directly or indirectly. Empty if this
            has not been precomputed.
        relation_func : function
            Function with arguments (node, graph) that generates objects
            with some relationship with node on the given graph.

        Returns
        -------
        bool
            True if t1 has the specified relationship with t2, either
            directly or through a series of intermediates; False otherwise.
        """
        # if id2 is None, or both are None, then it's by definition isa:
        if id2 is None or (id2 is None and id1 is None):
            return True
        # If only id1 is None, then it cannot be isa
        elif id1 is None:
            return False

        # If both terms are in the closure set then we can just look them
        # up and return
        if closure_set:
            term1 = self.get_uri(ns1, id1)
            term2 = self.get_uri(ns2, id2)
            if self._term_in_closure_namespace(term1) and \
                    self._term_in_closure_namespace(term2):
                return (term1, term2) in closure_set

        # Otherwise we do an actual graph query in the RDF graph
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
        return t2 in self.graph.transitiveClosure(relation_func, t1)

    def isa(self, ns1, id1, ns2, id2):
        """Return True if one entity has an "isa" relationship to another.

        Parameters
        ----------
        ns1 : str
            Namespace code for an entity.
        id1 : string
            URI for an entity.
        ns2 : str
            Namespace code for an entity.
        id2 : str
            URI for an entity.

        Returns
        -------
        bool
            True if t1 has an "isa" relationship with t2, either directly or
            through a series of intermediates; False otherwise.
        """
        rel_fun = lambda node, graph: self.isa_objects(node)
        return self.directly_or_indirectly_related(ns1, id1, ns2, id2,
                                                   self.isa_closure,
                                                   rel_fun)

    def partof(self, ns1, id1, ns2, id2):
        """Return True if one entity is "partof" another.

        Parameters
        ----------
        ns1 : str
            Namespace code for an entity.
        id1 : str
            URI for an entity.
        ns2 : str
            Namespace code for an entity.
        id2 : str
            URI for an entity.

        Returns
        -------
        bool
            True if t1 has a "partof" relationship with t2, either directly or
            through a series of intermediates; False otherwise.
        """
        rel_fun = lambda node, graph: self.partof_objects(node)
        return self.directly_or_indirectly_related(ns1, id1, ns2, id2,
                                                   self.partof_closure,
                                                   rel_fun)

    def isa_or_partof(self, ns1, id1, ns2, id2):
        """Return True if two entities are in an "isa" or "partof" relationship

        Parameters
        ----------
        ns1 : str
            Namespace code for an entity.
        id1 : str
            URI for an entity.
        ns2 : str
            Namespace code for an entity.
        id2 : str
            URI for an entity.

        Returns
        -------
        bool
            True if t1 has a "isa" or "partof" relationship with t2, either
            directly or through a series of intermediates; False otherwise.
        """
        rel_fun = lambda node, graph: self.isa_or_partof_objects(node)
        return self.directly_or_indirectly_related(ns1, id1, ns2, id2,
                                                   self.isa_or_partof_closure,
                                                   rel_fun)

    def get_equals(self, ns1, id1):
        u1 = self.get_uri(ns1, id1)
        t1 = rdflib.term.URIRef(u1)
        rel = rdflib.term.URIRef(self.relations_prefix + 'is_equal')
        to = [t.toPython() for t in list(self.graph.objects(t1, rel))]
        return to

    def get_opposites(self, ns1, id1):
        u1 = self.get_uri(ns1, id1)
        t1 = rdflib.term.URIRef(u1)
        rel = rdflib.term.URIRef(self.relations_prefix + 'is_opposite')
        to = [t.toPython() for t in list(self.graph.objects(t1, rel))]
        return to

    def is_opposite(self, ns1, id1, ns2, id2):
        """Return True if two entities are in an "is_opposite" relationship

        Parameters
        ----------
        ns1 : str
            Namespace code for an entity.
        id1 : str
            URI for an entity.
        ns2 : str
            Namespace code for an entity.
        id2 : str
            URI for an entity.

        Returns
        -------
        bool
            True if t1 has an "is_opposite" relationship with t2.
        """
        u2 = self.get_uri(ns2, id2)

        if u2 in self.get_opposites(ns1, id1):
            return True
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
        # First do a search in the set to see if there are any parents
        all_parents = {p for c, p in self.isa_or_partof_closure
                       if c == uri}
        # If there are no parents or we are looking for all, we can return here
        if not all_parents or type == 'all':
            return all_parents

        # If we need immediate parents, we search again, this time knowing that
        # the uri is definitely in the graph since it has some parents
        if type == 'immediate':
            node = rdflib.term.URIRef(uri)
            immediate_parents = list(set(self.isa_or_partof_objects(node)))
            return [p.toPython() for p in immediate_parents]
        elif type == 'top':
            # Here we iterate over all parents and find ones that have no
            # parents in the closure
            top_parents = [p for p in all_parents if
                           not {pp for pp, _ in self.isa_or_partof_closure
                                if pp == p}]
            return top_parents

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
            return 'http://identifiers.org/hgnc/' + id
        elif ns == 'UP':
            return 'http://identifiers.org/uniprot/' + id
        elif ns == 'FPLX':
            return 'http://identifiers.org/fplx/' + id
        elif ns == 'CHEBI':
            return 'http://identifiers.org/chebi/' + id
        elif ns in ['UN', 'WDI', 'FAO', 'HUME']:
            return \
                'https://github.com/clulab/eidos/wiki/JSON-LD/Grounding#' + id
        elif ns == 'SOFIA':
            return \
                'http://cs.cmu.edu/sofia/' + id
        elif ns == 'CWMS':
            if id.lower().startswith('ont::'):
                id = id[5:]
            return 'http://trips.ihmc.us/concepts/' + id.lower()
        elif ns == 'INDRA_ACTIVITIES':
            return 'http://sorger.med.harvard.edu/indra/activities/' + id
        elif ns == 'INDRA_MODS':
            return 'http://sorger.med.harvard.edu/indra/modifications/' + id
        elif ns == 'INDRA_LOCATIONS':
            return 'http://sorger.med.harvard.edu/indra/locations/' + id
        else:
            return ns + id

    @staticmethod
    def ns_id_from_uri(uri):
        sep_ix = uri.rfind('/') + 1
        ag_ns = uri[0:sep_ix]
        ag_id = uri[sep_ix:]
        ag_ns_name = ns_map.get(ag_ns)
        if ag_ns_name is None:
            raise UnknownNamespaceException('Unknown namespace %s' % ag_ns)
        return ag_ns_name, ag_id


class YamlHierarchyManager(HierarchyManager):
    """Class to manage YAML-based hierarchies.

    Parameters
    ----------
    root : dict
        A YAML data structure loaded with the yaml package.
    yaml_to_rdf : function
        A function that takes the root object as an argument and returns
        an RDF graph.
    """
    def __init__(self, root, yaml_to_rdf, add_leaves):
        self.yaml_root = root
        self.yaml_to_rdf = yaml_to_rdf
        self.add_leaves = add_leaves
        super(YamlHierarchyManager, self).__init__(None, True, True)
        G = self.yaml_to_rdf(self.yaml_root, self.add_leaves)
        self.load_from_rdf_graph(G)

    def add_entry(self, entry, examples=None):
        """Add a given entry to the ontology.

        Parameters
        ----------
        entry : str
            An entry in the ontology, with parts separated by /, e.g.,
            animals/mammals/dog.
        examples : Optional[list[str]]
            A list of strings that serve as examples for the given entry.
        """
        # TODO: Add the entry by finding the right place in the YAML object
        examples = examples if examples else []
        parts = entry.split('/')
        root = self.yaml_root
        for idx, part in enumerate(parts):
            new_root = None
            for element in root:
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

        G = self.yaml_to_rdf(self.yaml_root, self.add_leaves)
        self.load_from_rdf_graph(G)


def get_bio_hierarchies(from_pickle=True):
    """Return default hierarchies for the Bio context.

    Parameters
    ----------
    from_pickle : Optional[bool[
        If True, hierarchies are loded from a pre-generated pickle file.
        Otherwise, they are regenerated from RDF files (slower).
        Default: True

    Returns
    -------
    dict[str, HierarchyManager]
        A dict of hierarchy managers for each type of hierarchy.
    """
    if from_pickle:
        import pickle
        hierarchy_file = os.path.dirname(os.path.abspath(__file__)) + \
            '/../resources/bio_hierarchies.pkl'
        with open(hierarchy_file, 'rb') as fh:
            hierarchies = pickle.load(fh)
        return hierarchies

    def resource_path(fname):
        return os.path.join(os.path.dirname(__file__), os.pardir, 'resources',
                            fname)

    # Default entity hierarchy loaded from the RDF file at
    # `resources/entity_hierarchy.rdf`.
    entity_hierarchy = HierarchyManager(resource_path('entity_hierarchy.rdf'),
                                        build_closure=[
                                            'http://identifiers.org/hgnc',
                                            'http://identifiers.org/uniprot',
                                            'http://identifiers.org/fplx'
                                            ],
                                        uri_as_name=True)
    # Default modification hierarchy loaded from the RDF file at
    # `resources/modification_hierarchy.rdf`.
    modification_hierarchy = \
        HierarchyManager(resource_path('modification_hierarchy.rdf'),
                         build_closure=True, uri_as_name=True)
    # Default activity hierarchy loaded from the RDF file at
    # `resources/activity_hierarchy.rdf`.
    activity_hierarchy = \
        HierarchyManager(resource_path('activity_hierarchy.rdf'),
                         build_closure=True, uri_as_name=True)
    # Default cellular_component hierarchy loaded from the RDF file at
    # `resources/cellular_component_hierarchy.rdf`.
    ccomp_hierarchy = \
        HierarchyManager(resource_path('cellular_component_hierarchy.rdf'),
                         build_closure=False, uri_as_name=False)

    hierarchies = {'entity': entity_hierarchy,
                   'modification': modification_hierarchy,
                   'activity': activity_hierarchy,
                   'cellular_component': ccomp_hierarchy}
    return hierarchies


hierarchies = get_bio_hierarchies()


def get_wm_hierarchies():
    """Return default hierarchy managers for the World Modeling context.

    Returns
    -------
    dict[str, HierarchyManager]
        A dict of hierarchy managers for each type of hierarchy, in this context
        only an `entity` hierarchy is provided in the dict.
    """
    this_dir = os.path.dirname(__file__)
    wm_ont = os.path.join(this_dir, os.pardir, 'resources', 'wm_ontology.rdf')
    trips_ont = os.path.join(this_dir, os.pardir, 'sources', 'cwms',
                             'trips_ontology.rdf')
    sofia_ont = os.path.join(this_dir, os.pardir, 'sources', 'sofia',
                             'sofia_ontology.rdf')
    hm = HierarchyManager(wm_ont, build_closure=False, uri_as_name=True)
    hm.extend_with(trips_ont)
    hm.extend_with(sofia_ont)
    wm_hierarchies = {'entity': hm}
    return wm_hierarchies


class UnknownNamespaceException(Exception):
    pass
