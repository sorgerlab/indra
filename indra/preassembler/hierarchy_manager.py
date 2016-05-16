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

    @functools32.lru_cache(maxsize=1000)
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

# Load the default entity and modification hierarchies
entity_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/entity_hierarchy.rdf')
mod_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/modification_hierarchy.rdf')
act_file_path = os.path.join(os.path.dirname(__file__),
                    '../resources/activity_hierarchy.rdf')
entity_hierarchy = HierarchyManager(entity_file_path)
"""Default entity hierarchy loaded from the RDF file at
`resources/entity_hierarchy.rdf`."""
modification_hierarchy = HierarchyManager(mod_file_path)
"""Default modification hierarchy loaded from the RDF file at
`resources/modification_hierarchy.rdf`."""
activity_hierarchy = HierarchyManager(act_file_path)
"""Default activity hierarchy loaded from the RDF file at
`resources/activity_hierarchy.rdf`."""

