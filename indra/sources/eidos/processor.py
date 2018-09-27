from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import json
import logging
import objectpath
from indra.statements import Influence, Concept, Evidence


logger = logging.getLogger('eidos')

class EidosProcessor(object):
    """This processor extracts INDRA Statements from Eidos JSON-LD output.

    Parameters
    ----------
    json_dict : dict
        A JSON dictionary containing the Eidos extractions in JSON-LD format.

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
        self.sentence_dict = {}
        self.entity_dict = {}

    def get_events(self):
        extractions = \
            self.tree.execute("$.extractions[(@.@type is 'Extraction')]")
        if not extractions:
            return
        # Listify for multiple reuse
        extractions = list(extractions)

        events = [e for e in extractions if 'DirectedRelation' in
                  e.get('labels', [])]

        # Build a dictionary of entities and sentences by ID for convenient
        # lookup
        entities = [e for e in extractions if 'Concept' in
                    e.get('labels', [])]
        self.entity_dict = {entity['@id']: entity for entity in entities}

        documents = \
            self.tree.execute("$.documents[(@.@type is 'Document')]")
        self.sentence_dict = {}
        for document in documents:
            sentences = document.get('sentences', [])
            self.sentence_dict = {sent['@id']: sent for sent in sentences}

        # The first state corresponds to increase/decrease
        def get_polarity(x):
            # x is either subj or obj
            if 'states' in x:
                if x['states'][0]['type'] == 'DEC':
                    return -1
                elif x['states'][0]['type'] == 'INC':
                    return 1
                else:
                    return None
            else:
                return None

        def get_adjectives(x):
            # x is either subj or obj
            if 'states' in x:
                if 'modifiers' in x['states'][0]:
                    return [mod['text'] for mod in
                            x['states'][0]['modifiers']]
                else:
                    return []
            else:
                return []

        def _get_grounding_tuples(grounding):
            if not grounding or "values" not in grounding:
                return None

            grounding_values = grounding.get("values", [])

            grounding_tuples = map(
                # For some versions of eidos, groundings may erroneously have
                # the /examples suffix; strip that off if present
                lambda t: (
                    (t[0][: -len("/examples")], t[1])
                    if t[0].endswith("/examples")
                    else (t[0][: -len("/description")], t[1])
                    if t[0].endswith("/description")
                    else t
                ),
                map(
                    lambda x: (
                        (x["ontologyConcept"][1:], x["value"])
                        if x["ontologyConcept"].startswith("/")
                        else (x["ontologyConcept"], x["value"])
                    ),
                    # get all the groundings that have non-zero score
                    filter(lambda x: x["value"] > 0, grounding_values),
                ),
            )

            return list(grounding_tuples)

        def _get_groundings(entity):
            """Return Eidos groundings as a list of tuples with scores."""
            refs = {'TEXT': entity['text']}
            groundings = entity.get('groundings', [])
            for g in groundings:
                values = _get_grounding_tuples(g)
                # Only add these groundings if there are actual values listed
                if values:
                    key = g['name'].upper()
                    refs[key] = values
            return refs

        def _make_concept(entity):
            """Return Concept from an Eidos entity."""
            # Use the canonical name as the name of the Concept
            name = entity['canonicalName']
            # Save raw text and Eidos scored groundings as db_refs
            db_refs = {'TEXT': entity['text']}
            groundings = _get_groundings(entity)
            db_refs.update(groundings)
            concept = Concept(name, db_refs=db_refs)
            return concept

        def find_arg(event, arg_type):
            args = event.get('arguments', {})
            obj_tag = [arg for arg in args if arg['type'] == arg_type]
            if obj_tag:
                obj_id = obj_tag[0]['value']['@id']
            else:
                obj_id = None
            return obj_id

        for event in events:
            if not 'Causal' in event['labels']:
                continue

            # For now, just take the first source and first destination.
            # Later, might deal with hypergraph representation.
            subj_id = find_arg(event, 'source')
            obj_id = find_arg(event, 'destination')
            if subj_id is None or obj_id is None:
                continue

            subj = self.entity_dict[subj_id]
            obj = self.entity_dict[obj_id]

            subj_delta = {'adjectives': get_adjectives(subj),
                          'polarity': get_polarity(subj)}
            obj_delta = {'adjectives': get_adjectives(obj),
                         'polarity': get_polarity(obj)}

            evidence = self._get_evidence(event)

            st = Influence(_make_concept(subj), _make_concept(obj),
                           subj_delta, obj_delta, evidence=evidence)

            self.statements.append(st)

    def _get_evidence(self, event):
        """Return the Evidence object for the INDRA Statment."""
        provenance = event.get('provenance')

        # First try looking up the full sentence through provenance
        text = None
        time_annot = {}
        if provenance:
            sentence_tag = provenance[0].get('sentence')
            if sentence_tag and '@id' in sentence_tag:
                sentence_id = sentence_tag['@id']
                sentence = self.sentence_dict.get(sentence_id)
                if sentence is not None:
                    text = self._sanitize(sentence['text'])
                # Get temporal constraints if available
                timexes = sentence.get('timexes', [])
                if timexes:
                    time_text = timexes[0].get('text')
                    constraint = timexes[0]['intervals'][0]
                    start = None if constraint['start'] == 'Undef' else \
                        constraint['start']
                    end = None if constraint['end'] == 'Undef' else \
                        constraint['end']
                    duration = constraint['duration']
                    time_annot = {'text': time_text, 'start': start,
                                  'end': end, 'duration': duration}

        # If that fails, we can still get the text of the event
        if text is None:
            text = self._sanitize(event.get('text'))

        annotations = {
                'found_by'   : event.get('rule'),
                'provenance' : provenance,
                }
        if time_annot:
            annotations['time'] = time_annot
        ev = Evidence(source_api='eidos', text=text, annotations=annotations)
        return [ev]

    @staticmethod
    def _sanitize(text):
        """Return sanitized Eidos text field for human readability."""
        d = {'-LRB-': '(', '-RRB-': ')'}
        return re.sub('|'.join(d.keys()), lambda m: d[m.group(0)], text)
