import os
import rdflib
import functools32

def get_term(ns, id):
    if ns == 'HGNC':
        return 'http://identifiers.org/hgnc.symbol/' + id
    elif ns == 'UP':
        return 'http://identifiers.org/uniprot/' + id
    elif ns == 'BE' or ns == 'INDRA':
        return 'http://sorger.med.harvard.edu/indra/entities/' + id
    else:
        raise ValueError('Unknown namespace %s' % ns)

class HierarchyManager(object):
    """Store hierarchical relationships between different types of entities.

    Used to store, e.g., entity hierarchies (proteins and protein families)
    and modification hierarchies (serine phosphorylation vs. phosphorylation).

    Parameters
    ----------
    rdf_file : string
        Path to the RDF file containing the hierarchy.

    Attributes
    ----------
    graph : instance of `rdflib.Graph`
        The RDF graph containing the hierarchy.
    """
    prefixes = """
        PREFIX rn: <http://sorger.med.harvard.edu/indra/relations/>
        """

    def __init__(self, rdf_file):
        """Initialize with the path to an RDF file"""
        self.graph = rdflib.Graph()
        self.graph.parse(rdf_file)
        self.isa_closure = {}
        self.partof_closure = {}

    def build_transitive_closures(self):
        """Build the transitive closures of the hierarchy.

        This method constructs dictionaries which contain terms in the
        hierarchy as keys and either all the "isa+" or "partof+" related terms
        as values.
        """
        for rel, tc_dict in (('isa', self.isa_closure),
                             ('partof', self.partof_closure)):
            qstr = self.prefixes + """
                SELECT ?x ?y WHERE {{
                    {{?x rn:{0}+ ?y .}}
                    }}
                """.format(rel)
            res = self.graph.query(qstr)
            for x, y in res:
                xs = x.toPython()
                ys = y.toPython()
                try:
                    tc_dict[xs].append(ys)
                except KeyError:
                    tc_dict[xs] = [ys]

    """
    @functools32.lru_cache(maxsize=100000)
    def find_entity(self, x):
        Get the entity that has the specified name (or synonym).

        Parameters
        ----------
        x : string
            Name or synonym for the target entity.

        qstr = self.prefixes +
            SELECT ?x WHERE {{
                {{
                {{ ?x rn:hasName "{0}" . }}
                UNION
                {{ ?x rn:hasSynonym "{0}" . }}
                UNION
                {{ ?x rn:hasId "{0}" . }}
                }}
            }}
            .format(x)
        res = self.graph.query(qstr)
        if list(res):
            en = list(res)[0][0].toPython()
            return en
        else:
            return None
    """

    @functools32.lru_cache(maxsize=100000)
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
        term1 = get_term(ns1, id1)
        term2 = get_term(ns2, id2)

        if self.isa_closure:
            ec = self.isa_closure.get(term1)
            if ec is not None and term2 in ec:
                return True
            else:
                return False

        qstr = self.prefixes + """ 
            SELECT (COUNT(*) as ?s) WHERE {{
                <{}> rn:isa+ <{}> .
                }}
            """.format(term1, term2)
        res = self.graph.query(qstr)
        count = [r[0] for r in res][0]
        if count.toPython() == 1:
            return True
        else:
            return False

    @functools32.lru_cache(maxsize=100000)
    def partof(self, t1, t2):
        """Indicate whether one entity is physically part of another.

        Parameters
        ----------
        t1 : string
            URI for an entity.
        t2 : string
            URI for an entity.

        Returns
        -------
        bool
            True if t1 has a "partof" relationship with t2, either directly or
            through a series of intermediates; False otherwise.
        """
        en1 = self.find_entity(t1)
        en2 = self.find_entity(t2)

        if en1 is None or en2 is None:
            return None

        if self.transitive_closure:
            ec = self.transitive_closure.get(en1)
            if ec and en2 in ec:
                return True
            else:
                return False

        qstr = self.prefixes + """ 
            SELECT (COUNT(*) as ?s) WHERE {{
                <{}> rn:partof+ <{}> .
                }}
            """.format(en1, en2)
        res = self.graph.query(qstr)
        count = [r[0] for r in res][0]
        if count.toPython() == 1:
            return True
        else:
            return False

# Load the default entity and modification hierarchies
entity_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/entity_hierarchy.rdf')
mod_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/modification_hierarchy.rdf')
act_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/activity_hierarchy.rdf')
ccomp_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/cellular_component_hierarchy.rdf')
"""Default entity hierarchy loaded from the RDF file at
`resources/entity_hierarchy.rdf`."""
entity_hierarchy = HierarchyManager(entity_file_path)
entity_hierarchy.build_transitive_closures()
"""Default modification hierarchy loaded from the RDF file at
`resources/modification_hierarchy.rdf`."""
modification_hierarchy = HierarchyManager(mod_file_path)
modification_hierarchy.build_transitive_closures()
"""Default activity hierarchy loaded from the RDF file at
`resources/activity_hierarchy.rdf`."""
activity_hierarchy = HierarchyManager(act_file_path)
activity_hierarchy.build_transitive_closures()
"""Default cellular_component hierarchy loaded from the RDF file at
`resources/cellular_component_hierarchy.rdf`."""
ccomp_hierarchy = HierarchyManager(ccomp_file_path)

hierarchies = {'entity': entity_hierarchy,
               'modification': modification_hierarchy,
               'activity': activity_hierarchy,
               'cellular_component': ccomp_hierarchy}
