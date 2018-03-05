import rdflib
import collections
from indra.statements import Agent, Influence
from rdflib.serializer import Serializer

prefixes = """
    PREFIX causal: <http://worldmodelers.com/CauseEffect#>
    PREFIX cco: <http://www.ontologyrepository.com/CommonCoreOntologies/>
    """

class BBNProcessor(object):
    """Process a BBN JSON-LD file into INDRA statements.

    Parameters
    ----------
    fname: str
        JSON-LD filename with BBN reader output

    Attributes
    ----------
    stmts: list<indra.statements.Statement>
        INDRA statements extracted from BBN reader output
    """
    def __init__(self, fname):
        self.stmts = self.process_json_file(fname)

    def load_graph(self, fname):
        g = rdflib.Graph()
        with open(fname, 'rb') as fh:
            print('Started loading graph')
            g.parse(fh, format='json-ld')
            print('Finished loading graph')
        return g

    def process_json_file(self, fname):
        g = self.load_graph(fname)
        query = prefixes + """
            SELECT ?rel
                ?causetext
                ?effecttext
                ?evtext
            WHERE {
                ?rel a causal:CausalAssertion .
                ?rel causal:has_cause ?cause .
                ?rel causal:has_effect ?effect .
                ?rel cco:has_text_value ?evtext .
                ?cause cco:has_text_value ?causetext .
                ?effect cco:has_text_value ?effecttext .
                }
            """
        # All this stuff below is just to get the shorter of the two text
        # values for the cause/effect events
        res = g.query(query)
        rdict = collections.defaultdict(list)
        for rel, causetext, effecttext, evtext in res:
            relid = rel.rsplit('#')[1]
            rdict[relid].append((causetext, effecttext))
        stmts = []
        for relid, ces in rdict.items():
            cause = sorted(set([c[0] for c in ces]), key=lambda x: len(x))[0]
            effect = sorted(set([c[1] for c in ces]), key=lambda x: len(x))[0]
            print('%s -> %s' % (cause, effect))
            stmt = Influence(Agent(cause), Agent(effect))
            stmts.append(stmt)
        return stmts


