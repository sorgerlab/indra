import rdflib
import logging
import collections
from indra.statements import Agent, Influence


logger = logging.getLogger('bbn')


prefixes = """
    PREFIX causal: <http://worldmodelers.com/CauseEffect#>
    PREFIX cco: <http://www.ontologyrepository.com/CommonCoreOntologies/>
    """


class BBNProcessor(object):
    """Process a BBN extraction graph into INDRA Statements.

    Parameters
    ----------
    graph : rdflib.Graph
        An rdflib graph representing extractions, typically loaded from
        a JSON-LD file.

    Attributes
    ----------
    statements: list[indra.statements.Statement]
        INDRA statements extracted from BBN reader output.
    """
    def __init__(self, graph):
        self.graph = graph
        self.statements = []

    def get_statements(self):
        # SPARQL query to get causal relations and their arguments
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
        # Run the query
        res = self.graph.query(query)
        # All this stuff below is just to get the shorter of the two text
        # values for the cause/effect events
        rdict = collections.defaultdict(list)
        for rel, causetext, effecttext, evtext in res:
            relid = rel.rsplit('#')[1]
            rdict[relid].append((causetext, effecttext))
        for relid, ces in rdict.items():
            cause = sorted(set([c[0] for c in ces]), key=lambda x: len(x))[0]
            effect = sorted(set([c[1] for c in ces]), key=lambda x: len(x))[0]
            logger.debug('%s -> %s' % (cause, effect))
            stmt = Influence(Agent(cause), Agent(effect))
            self.statements.append(stmt)
