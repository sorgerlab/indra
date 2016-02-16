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

    def find_entity(self, x):
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
        
    def isa(self, t1, t2):
        """Both t1 and t2 are entities and we determine whether there is a
        series of "isa" edges from t1 to t2.
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
