import rdflib

class HierarchyManager(object):
    prefixes = """
        PREFIX rn: <http://sorger.med.harvard.edu/indra/relations/>
        PREFIX en: <http://sorger.med.harvard.edu/indra/entities/>
        """

    def __init__(self, rdf_file):
        """Initialize with the path to an RDF file"""
        self.graph = rdflib.Graph()
        self.graph.parse(rdf_file)

    def isa(self, t1, t2):
        """Both t1 and t2 are entities and we determine whether there is a
        series of "isa" edges from t1 to t2.
        """
        qstr = self.prefixes + """ 
            SELECT (COUNT(*) as ?s) WHERE {{
                en:{} rn:isa+ en:{} .
                }}
            """.format(t1, t2)
        res = self.graph.query(qstr)
        count = [r[0] for r in res][0]
        if count.toPython() == 1:
            return True
        else:
            return False
