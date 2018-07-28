import re
import os
import rdflib
import logging
import objectpath
import collections
from os.path import basename
from indra.statements import Concept, Influence, Evidence


logger = logging.getLogger('hume')


class HumeJsonLdProcessor(object):
    """This processor extracts INDRA Statements from Hume JSON-LD output.

    Parameters
    ----------
    json_dict : dict
        A JSON dictionary containing the Hume extractions in JSON-LD format.

    Attributes
    ----------
    tree : objectpath.Tree
        The objectpath Tree object representing the extractions.
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted by the processor.
    """
    def __init__(self, json_dict):
        self.tree = objectpath.Tree(json_dict)
        self.statements = []
        self.document_dict = {}
        self.concept_dict = {}
        self.relation_dict = {}
        self.eid_stmt_dict = {}

    def get_documents(self):
        """Populate the sentences attribute with a dict keyed by document id."""
        documents = self.tree.execute("$.documents")
        for doc in documents:
            sentences = {s['@id']: s['text'] for s in doc.get('sentences', [])}
            self.document_dict[doc['@id']] = {'sentences': sentences,
                                              'location': doc['location']}
        return

    def get_events(self):
        # Get all extractions
        extractions = \
            list(self.tree.execute("$.extractions[(@.@type is 'Extraction')]"))

        # Get relations from extractions
        relations = [e for e in extractions if 'DirectedRelation' in
                     e.get('labels', [])]
        if not relations:
            return
        self.relation_dict = {rel['@id']: rel for rel in relations}

        # List out relation types and their default (implied) polarities.
        relation_polarities = {'causation': 1, 'precondition': 1, 'catalyst': 1,
                               'mitigation': -1, 'prevention': -1}

        # Restrict to known relation types
        relations = [r for r in relations if any([rt in r.get('type') for rt in
                                                  relation_polarities.keys()])]
        logger.info('%d relations of types %s found'
                    % (len(relations), ', '.join(relation_polarities.keys())))

        # Build a dictionary of concepts and sentences by ID for convenient
        # lookup
        concepts = [e for e in extractions if
                    set(e.get('labels', [])) & {'Event', 'Entity'}]
        self.concept_dict = {concept['@id']: concept for concept in concepts}

        self.get_documents()

        for relation in relations:
            relation_type = relation.get('type')
            subj_concept, subj_delta = self._get_concept(relation, 'source')
            obj_concept, obj_delta = self._get_concept(relation, 'destination')

            # Apply the naive polarity from the type of statement. For the
            # purpose of the multiplication here, if obj_delta['polarity'] is
            # None to begin with, we assume it is positive
            obj_delta['polarity'] = \
                relation_polarities[relation_type] * \
                (obj_delta['polarity'] if obj_delta['polarity'] is not None
                 else 1)

            if not subj_concept or not obj_concept:
                continue

            evidence = self._get_evidence(relation, subj_concept, obj_concept,
                                          get_states(relation))

            st = Influence(subj_concept, obj_concept, subj_delta, obj_delta,
                           evidence=evidence)
            self.eid_stmt_dict[relation['@id']] = st
            self.statements.append(st)

        return

    def _make_concept(self, entity):
        """Return Concept from a Hume entity."""
        # Use the canonical name as the name of the Concept by default
        name = self._sanitize(entity['canonicalName'])
        # But if there is a trigger head text, we prefer that since
        # it almost always results in a cleaner name
        # This is removed for now since the head word seems to be too
        # minimal for some concepts, e.g. it gives us only "security"
        # for "food security".
        """
        trigger = entity.get('trigger')
        if trigger is not None:
            head_text = trigger.get('head text')
            if head_text is not None:
                name = head_text
        """
        # Save raw text and Hume scored groundings as db_refs
        db_refs = {'TEXT': entity['text']}
        hume_grounding = _get_hume_grounding(entity)
        # We could get an empty list here in which case we don't add the
        # grounding
        if hume_grounding:
            db_refs['HUME'] = hume_grounding
        concept = Concept(name, db_refs=db_refs)
        return concept

    def _get_concept(self, event, arg_type):
        eid = _choose_id(event, arg_type)
        ev = self.concept_dict[eid]
        concept = self._make_concept(ev)
        ev_delta = {'adjectives': [],
                    'states': get_states(ev),
                    'polarity': get_polarity(ev)}
        return concept, ev_delta

    def _get_evidence(self, event, subj_concept, obj_concept, adjectives):
        """Return the Evidence object for the INDRA Statement."""
        provenance = event.get('provenance')

        # First try looking up the full sentence through provenance
        doc_id = provenance[0]['document']['@id']
        sent_id = provenance[0]['sentence']
        text = self.document_dict[doc_id]['sentences'][sent_id]
        text = self._sanitize(text)
        bounds = [provenance[0]['documentCharPositions'][k]
                  for k in ['start', 'end']]

        annotations = {
            'found_by': event.get('rule'),
            'provenance': provenance,
            'event_type': basename(event.get('type')),
            'adjectives': adjectives,
            'bounds': bounds
            }
        location = self.document_dict[doc_id]['location']
        ev = Evidence(source_api='hume', text=text, annotations=annotations,
                      pmid=location)
        return [ev]

    @staticmethod
    def _sanitize(text):
        """Return sanitized Hume text field for human readability."""
        # TODO: any cleanup needed here?
        if text is None:
            return None
        text = text.replace('\n', ' ')
        return text


def _choose_id(event, arg_type):
    args = event.get('arguments', {})
    obj_tag = [arg for arg in args if arg['type'] == arg_type]
    if obj_tag:
        obj_id = obj_tag[0]['value']['@id']
    else:
        obj_id = None
    return obj_id


def get_states(event):
    ret_list = []
    if 'states' in event:
        for state_property in event['states']:
            if state_property['type'] != 'polarity':
                ret_list.append(state_property['text'])
    return ret_list


def _get_hume_grounding(entity):
    """Return Hume grounding."""
    groundings = entity.get('grounding')
    if not groundings:
        return None
    def get_ont_concept(concept):
        """Strip slah, replace spaces and remove example leafs."""
        if concept.startswith('/'):
            concept = concept[1:]
        concept = concept.replace(' ', '_')
        # We eliminate any entries that aren't ontology categories
        # these are typically "examples" corresponding to the category
        while concept not in hume_onto_entries:
            parts = concept.split('/')
            if len(parts) == 1:
                break
            concept = '/'.join(parts[:-1])
        return concept

    # Basic collection of grounding entries
    raw_grounding_entries = [(get_ont_concept(g['ontologyConcept']),
                              g['value']) for g in groundings]

    # Occasionally we get duplicate grounding entries, we want to
    # eliminate those here
    grounding_dict = {}
    for cat, score in raw_grounding_entries:
        if (cat not in grounding_dict) or (score > grounding_dict[cat]):
            grounding_dict[cat] = score
    # Then we sort the list in reverse order according to score
    # Sometimes the exact same score appears multiple times, in this
    # case we prioritize by the "depth" of the grounding which is
    # obtained by looking at the number of /-s in the entry
    grounding_entries = sorted(list(set(grounding_dict.items())),
                               key=lambda x: (x[1], x[0].count('/')),
                               reverse=True)

    return grounding_entries


def get_polarity(event):
    pol_map = {'Positive': 1, 'Negative': -1}
    if 'states' in event:
        for state_property in event['states']:
            if state_property['type'] == 'polarity':
                return pol_map[state_property['text']]
    return None


def _get_ontology_entries():
    path_here = os.path.dirname(os.path.abspath(__file__))
    onto_file = os.path.join(path_here, 'hume_ontology.rdf')
    G = rdflib.Graph()
    G.parse(onto_file, format='nt')
    entries = [e.toPython().split('#')[1] for e in G.all_nodes()]
    return entries


hume_onto_entries = _get_ontology_entries()


# OLD PROCESSOR

prefixes = """
    PREFIX causal: <http://www.bbn.com/worldmodelers/ontology/wm/CauseEffect#>
    PREFIX ev: <http://www.bbn.com/worldmodelers/ontology/wm/Event#>
    PREFIX prov: <http://www.bbn.com/worldmodelers/ontology/wm/DataProvenance#>
    PREFIX cco: <http://www.ontologyrepository.com/CommonCoreOntologies/>
    """


class HumeProcessor(object):
    """Process a Hume extraction graph into INDRA Statements.

    Parameters
    ----------
    graph : rdflib.Graph
        An rdflib graph representing extractions, typically loaded from
        a JSON-LD file.

    Attributes
    ----------
    statements: list[indra.statements.Statement]
        INDRA statements extracted from Hume reader output.
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
    is extracted from the Hume JSON file.

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
        ev = Evidence(source_api='hume', text=str(evidence_text))

        # Convert from rdf literal to python string
        cause_text = str(cause_text)
        effect_text = str(effect_text)

        # Make cause concept
        cause_db_refs = {'TEXT': cause_text}
        if self.cause_type is not None:
            cause_db_refs['HUME'] = self.cause_type
        cause_concept = Concept(cause_text, db_refs=cause_db_refs)

        # Make effect concept
        effect_db_refs = {'TEXT': effect_text}
        if self.effect_type is not None:
            effect_db_refs['HUME'] = self.effect_type
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
