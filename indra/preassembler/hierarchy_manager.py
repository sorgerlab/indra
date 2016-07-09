import os
import rdflib
import functools32

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
        PREFIX en: <http://sorger.med.harvard.edu/indra/entities/>
        """

    def __init__(self, rdf_file):
        """Initialize with the path to an RDF file"""
        self.graph = rdflib.Graph()
        self.graph.parse(rdf_file)
        self.transitive_closure = {}

    def build_transitive_closure(self, rel='isa'):
        """Build the transitive closure of the hierarchy.

        This method constructs a dictionary which contains terms in the
        hierarchy as keys and all the "isa+" related terms as values.

        Parameters
        ----------
        rel : string
            The relationship along which to build the transitive closure.
            By default this is "isa" but it can also be "partof".
        """
        qstr = self.prefixes + """
            SELECT ?x ?y WHERE {{
                {{?x rn:{0}+ ?y .}}
                }}
            """.format(rel)
        res = self.graph.query(qstr)
        self.transitive_closure = {}
        for x, y in res:
            xs = x.toPython()
            ys = y.toPython()
            try:
                self.transitive_closure[xs].append(ys)
            except KeyError:
                self.transitive_closure[xs] = [ys]

    @functools32.lru_cache(maxsize=100000)
    def find_entity(self, x):
        """Get the entity that has the specified name (or synonym).

        Parameters
        ----------
        x : string
            Name or synonym for the target entity.
        """
        qstr = self.prefixes + """
            SELECT ?x WHERE {{
                {{
                {{ ?x rn:hasName "{0}" . }}
                UNION
                {{ ?x rn:hasSynonym "{0}" . }}
                UNION
                {{ ?x rn:hasId "{0}" . }}
                }}
            }}
            """.format(x)
        res = self.graph.query(qstr)
        if list(res):
            en = list(res)[0][0].toPython()
            return en
        else:
            return None

    @functools32.lru_cache(maxsize=100000)
    def isa(self, t1, t2):
        """Indicate whether one entity has an "isa" relationship to another.

        Parameters
        ----------
        t1 : string
            URI for an entity.
        t2 : string
            URI for an entity.

        Returns
        -------
        bool
            True if t1 has an "isa" relationship with t2, either directly or
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
                <{}> rn:isa+ <{}> .
                }}
            """.format(en1, en2)
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
entity_hierarchy.build_transitive_closure()
"""Default modification hierarchy loaded from the RDF file at
`resources/modification_hierarchy.rdf`."""
modification_hierarchy = HierarchyManager(mod_file_path)
modification_hierarchy.build_transitive_closure()
"""Default activity hierarchy loaded from the RDF file at
`resources/activity_hierarchy.rdf`."""
activity_hierarchy = HierarchyManager(act_file_path)
activity_hierarchy.build_transitive_closure()
"""Default cellular_component hierarchy loaded from the RDF file at
`resources/cellular_component_hierarchy.rdf`."""
ccomp_hierarchy = HierarchyManager(ccomp_file_path)

hierarchies = {'entity': entity_hierarchy,
               'modification': modification_hierarchy,
               'activity': activity_hierarchy,
               'cellular_component': ccomp_hierarchy}
