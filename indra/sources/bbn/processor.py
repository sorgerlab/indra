import rdflib
import logging
import collections
from indra.statements import Concept, Influence, Evidence


logger = logging.getLogger('bbn')


prefixes = """
    PREFIX causal: <http://www.bbn.com/worldmodelers/ontology/wm/CauseEffect#>
    PREFIX ev: <http://www.bbn.com/worldmodelers/ontology/wm/Event#>
    PREFIX prov: <http://www.bbn.com/worldmodelers/ontology/wm/DataProvenance#>
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
        """Extract causal assertions in the graph into INDRA statements."""
        # SPARQL query to get causal relations and their arguments
        query = prefixes + """
            SELECT ?rel
                ?causetext
                ?effecttext
                ?evtext
                ?cause_polarity
                ?effect_polarity
                ?cause_type
                ?effect_type
            WHERE {
                ?rel a causal:CausalAssertion .
                ?rel causal:has_cause ?cause .
                ?rel causal:has_effect ?effect .
                ?rel prov:has_text_value ?evtext .
                ?cause prov:has_text_value ?causetext .
                ?effect prov:has_text_value ?effecttext .
                OPTIONAL
                    {?cause ev:has_polarity ?cause_polarity .
                    ?effect ev:has_polarity ?effect_polarity .
                    ?cause a ?cause_type .
                    ?effect a ?effect_type .}
                }
            """
        # Run the query
        res = self.graph.query(query)

        # Accumulate the cause, effect, and evidence textss for each causal
        # assertion. When there are several cause texts, the CauseEffect
        # class will only keep the shortest when we generate the statement
        # (since these events include both the full event and the cause/effect
        # snippet).
        rdict = collections.defaultdict(CauseEffect)
        for rel, cause_text, effect_text, evtext, cause_polarity, \
                effect_polarity, cause_type, effect_type in res:
            relid = shorter_name(rel)

            rdict[relid].cause_texts.add(cause_text)
            rdict[relid].effect_texts.add(effect_text)
            rdict[relid].evidence_texts.add(evtext)

            if cause_polarity is not None:
                rdict[relid].cause_polarity = shorter_name(cause_polarity)
            if effect_polarity is not None:
                rdict[relid].effect_polarity = shorter_name(effect_polarity)
            if cause_type is not None:
                rdict[relid].cause_type = shorter_name(cause_type)
            if effect_type is not None:
                rdict[relid].effect_type = shorter_name(effect_type)
        not_positive = 0
        for relid, ces in rdict.items():
            statement = ces.to_statement()
            if statement is None:  # Returns None when polarity not positive
                not_positive = not_positive + 1
            else:
                self.statements.append(statement)
        if not_positive > 0:
            print('%d statements skipped because of polarity' % not_positive)


class CauseEffect(object):
    """A data structure to incrementally store cause/effect information as it
    is extracted from the BBN JSON file.

    Parameters
    ----------
    cause_texts: list<str>
        A list of causes, in text
    effect_texts: list<str>
        A list of effects, in text
    evidence_texts: list<str>
        A list of evidence texts
    cause_polarity: str
        Polarity of the cause (no statement generated if not Positive)
    effect_polarity: str
        Polarity of the effect (no statement generated if not Positive)
    cause_type: str
        The type of cause
    effect_type: str
        The type of effect
    """

    def __init__(self):
        """Initialize with all fields blank; fields are populated as they
        are read from the JSON."""
        self.cause_texts = set()
        self.effect_texts = set()
        self.evidence_texts = set()
        self.cause_polarity = None
        self.effect_polarity = None
        self.cause_type = None
        self.effect_type = None

    def __repr__(self):
        """Convert to string representation, suitable for debugging."""
        ct = shortest_string_in_list(self.cause_texts)
        et = shortest_string_in_list(self.effect_texts)
        ev = ','.join(self.evidence_texts)
        return '%s -> %s [%s, %s, %s]' % (ct, et, ev,
                                          repr(self.cause_polarity),
                                          repr(self.effect_polarity))

    def to_statement(self):
        """Converts to an INDRA statement, or returns None if either the cause
        polarity or effect polarity is not Positive."""
        if self.cause_polarity != 'Positive' or \
                self.effect_polarity != 'Positive':
                    return None

        # The cause and effect events list both the full text and the text
        # identified as the cause/effect. Get the relevant text by getting
        # the shortest string.
        cause_text = shortest_string_in_list(self.cause_texts)
        effect_text = shortest_string_in_list(self.effect_texts)

        # Add an evidence object with the full text. There should be exactly
        # only full text string, but if there is more than one, list them all.
        # Note how we're careful to convert from rdflib's string representation
        # to a python string with str().
        evidence_texts = list(self.evidence_texts)
        if len(evidence_texts) == 1:
            evidence_text = evidence_texts[0]
        else:
            evidence_text = repr(evidence_texts)
        ev = Evidence(source_api='bbn', text=str(evidence_text))

        # Convert from rdf literal to python string
        cause_text = str(cause_text)
        effect_text = str(effect_text)

        # Make cause concept
        cause_db_refs = {'TEXT': cause_text}
        if self.cause_type is not None:
            cause_db_refs['BBN'] = self.cause_type
        cause_concept = Concept(cause_text, db_refs=cause_db_refs)

        # Make effect concept
        effect_db_refs = {'TEXT': effect_text}
        if self.effect_type is not None:
            effect_db_refs['BBN'] = self.effect_type
        effect_concept = Concept(effect_text, db_refs=effect_db_refs)

        return Influence(cause_concept, effect_concept, evidence=ev)


def shortest_string_in_list(string_list):
    """Given a list of strings, returns the shortest."""
    shortest_length = None
    shortest_string = None

    for s in string_list:
        if shortest_string is None or len(s) < shortest_length:
            shortest_length = len(s)
            shortest_string = s
    return shortest_string


def shorter_name(key):
    """Finds a shorter name for an id by only taking the last part of the URI,
    after the last / and the last #. Also replaces - and . with _.

    Parameters
    ----------
    key: str
        Some URI

    Returns
    -------
    key_short: str
        A shortened, but more ambiguous, identifier
    """
    key_short = key
    for sep in ['#', '/']:
        ind = key_short.rfind(sep)
        if ind is not None:
            key_short = key_short[ind+1:]
        else:
            key_short = key_short
    return key_short.replace('-', '_').replace('.', '_')
