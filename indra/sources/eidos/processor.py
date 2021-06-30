import re
import copy
import logging
import objectpath
from indra.statements import *


logger = logging.getLogger(__name__)


class EidosProcessor(object):
    """This processor extracts INDRA Statements from Eidos JSON-LD output.

    Parameters
    ----------
    json_dict : dict
        A JSON dictionary containing the Eidos extractions in JSON-LD format.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted by the processor.
    """
    def __init__(self, json_dict):
        self.doc = EidosDocument(json_dict)
        self.statements = []

    def extract_causal_relations(self):
        """Extract causal relations as Statements."""
        # Get the extractions that are labeled as directed and causal
        relations = [e for e in self.doc.extractions if
                     'DirectedRelation' in e['labels'] and
                     'Causal' in e['labels']]
        # For each relation, we try to extract an INDRA Statement and
        # save it if its valid
        for relation in relations:
            stmt = self.get_causal_relation(relation)
            if stmt is not None:
                self.statements.append(stmt)

    def extract_correlations(self):
        events = [e for e in self.doc.extractions if
                  'UndirectedRelation' in e['labels'] and
                  'Correlation' in e['labels']]
        for event in events:
            # For now, just take the first source and first destination.
            # Later, might deal with hypergraph representation.
            arg_ids = find_args(event, 'argument')
            if len(arg_ids) != 2:
                logger.warning('Skipping correlation with not 2 arguments.')

            # Resolve coreferences by ID
            arg_ids = [self.doc.coreferences.get(arg_id, arg_id)
                       for arg_id in arg_ids]

            # Get the actual entities
            args = [self.doc.entities[arg_id] for arg_id in arg_ids]
            # Make Events from the entities
            members = [self.get_event(arg) for arg in args]
            # Get the evidence
            evidence = self.get_evidence(event)

            st = Association(members, evidence=[evidence])
            self.statements.append(st)

    def extract_events(self):
        events = [e for e in self.doc.extractions if
                  'Concept-Expanded' in e['labels']]
        for event_entry in events:
            event = self.get_event(event_entry)
            evidence = self.get_evidence(event_entry)
            event.evidence = [evidence]
            if not event.context and evidence.context:
                event.context = copy.deepcopy(evidence.context)
                evidence.context = None
            self.statements.append(event)

    def extract_all_events(self):
        events = [e for e in self.doc.extractions if
                  {'Concept', 'Concept-Expanded'} & set(e['labels'])]
        for event_entry in events:
            event = self.get_event(event_entry)
            evidence = self.get_evidence(event_entry)
            event.evidence = [evidence]
            if not event.context and evidence.context:
                event.context = copy.deepcopy(evidence.context)
                evidence.context = None
            self.statements.append(event)


    def get_event_by_id(self, event_id):
        # Resolve coreferences by ID
        event_id = self.doc.coreferences.get(event_id, event_id)
        # Get the actual entity
        event = self.doc.entities[event_id]
        return self.get_event(event)

    def get_event(self, event):
        concept = self.get_concept(event)
        states = event.get('states', [])
        extracted_states = self.extract_entity_states(states)
        polarity = extracted_states.get('polarity')
        adjectives = extracted_states.get('adjectives')
        delta = QualitativeDelta(polarity=polarity, adjectives=adjectives)
        stmt = Event(concept, delta=delta)
        return stmt

    def get_causal_relation(self, relation):
        # For now, just take the first source and first destination.
        # Later, might deal with hypergraph representation.
        subj_id = find_arg(relation, 'source')
        obj_id = find_arg(relation, 'destination')
        if subj_id is None or obj_id is None:
            return None

        subj = self.get_event_by_id(subj_id)
        obj = self.get_event_by_id(obj_id)

        evidence = self.get_evidence(relation)

        # We also put the adjectives and polarities into annotations since
        # they could otherwise get squashed upon preassembly
        evidence.annotations['subj_polarity'] = subj.delta.polarity
        evidence.annotations['obj_polarity'] = obj.delta.polarity
        evidence.annotations['subj_adjectives'] = subj.delta.adjectives
        evidence.annotations['obj_adjectives'] = obj.delta.adjectives
        evidence.annotations['subj_context'] = subj.context.to_json() if \
            subj.context else {}
        evidence.annotations['obj_context'] = obj.context.to_json() if \
            obj.context else {}

        st = Influence(subj, obj, evidence=[evidence])
        return st

    def get_evidence(self, relation):
        """Return the Evidence object for the INDRA Statment."""
        provenance = relation.get('provenance')

        # First try looking up the full sentence through provenance
        text = None
        context = None
        if provenance:
            sentence_tag = provenance[0].get('sentence')
            if sentence_tag and '@id' in sentence_tag:
                sentence_id = sentence_tag['@id']
                sentence = self.doc.sentences.get(sentence_id)
                if sentence is not None:
                    text = _sanitize(sentence['text'])

            # Here we try to get the title of the document and set it
            # in the provenance
            doc_id = provenance[0].get('document', {}).get('@id')
            if doc_id:
                title = self.doc.documents.get(doc_id, {}).get('title')
                if title:
                    provenance[0]['document']['title'] = title

        annotations = {'found_by': relation.get('rule'),
                       'provenance': provenance}
        if self.doc.dct is not None:
            annotations['document_creation_time'] = self.doc.dct.to_json()

        epistemics = {}
        negations = self.get_negation(relation)
        hedgings = self.get_hedging(relation)
        if hedgings:
            epistemics['hedgings'] = hedgings
        if negations:
            # This is the INDRA standard to show negation
            epistemics['negated'] = True
            # But we can also save the texts associated with the negation
            # under annotations, just in case it's needed
            annotations['negated_texts'] = negations

        # If that fails, we can still get the text of the relation
        if text is None:
            text = _sanitize(relation.get('text'))

        ev = Evidence(source_api='eidos', text=text, annotations=annotations,
                      context=context, epistemics=epistemics)
        return ev

    @staticmethod
    def get_negation(event):
        """Return negation attached to an event.

        Example: "states": [{"@type": "State", "type": "NEGATION",
                             "text": "n't"}]
        """
        states = event.get('states', [])
        if not states:
            return []
        negs = [state for state in states
                if state.get('type') == 'NEGATION']
        neg_texts = [neg['text'] for neg in negs]
        return neg_texts

    @staticmethod
    def get_hedging(event):
        """Return hedging markers attached to an event.

        Example: "states": [{"@type": "State", "type": "HEDGE",
                             "text": "could"}
        """
        states = event.get('states', [])
        if not states:
            return []
        hedgings = [state for state in states
                    if state.get('type') == 'HEDGE']
        hedging_texts = [hedging['text'] for hedging in hedgings]
        return hedging_texts

    def extract_entity_states(self, states):
        if states is None:
            return {'polarity': None, 'adjectives': []}
        polarity = None
        adjectives = []
        for state in states:
            if polarity is None:
                if state['type'] == 'DEC':
                    polarity = -1
                    # Handle None entry here
                    mods = state.get('modifiers') if \
                        state.get('modifiers') else []
                    adjectives += [mod['text'] for mod in mods]
                elif state['type'] == 'INC':
                    polarity = 1
                    mods = state.get('modifiers') if \
                        state.get('modifiers') else []
                    adjectives += [mod['text'] for mod in mods]
                elif state['type'] == 'QUANT':
                    adjectives.append(state['text'])
        return {'polarity': polarity, 'adjectives': adjectives}

    def get_groundings(self, entity):
        """Return groundings as db_refs for an entity."""
        return {'TEXT': entity['text']}

    def get_concept(self, entity):
        """Return Concept from an Eidos entity."""
        # Use the canonical name as the name of the Concept
        name = entity['canonicalName']
        db_refs = self.get_groundings(entity)
        concept = Concept(name, db_refs=db_refs)
        return concept

    def get_all_events(self):
        """Return a list of all standalone events from a list
        of statements."""
        events = []
        for stmt in self.statements:
            stmt = copy.deepcopy(stmt)
            if isinstance(stmt, Influence):
                for member in [stmt.subj, stmt.obj]:
                    member.evidence = stmt.evidence[:]
                    # Remove the context since it may be for the other member
                    for ev in member.evidence:
                        ev.context = None
                    events.append(member)
            elif isinstance(stmt, Association):
                for member in stmt.members:
                    member.evidence = stmt.evidence[:]
                    # Remove the context since it may be for the other member
                    for ev in member.evidence:
                        ev.context = None
                    events.append(member)
            elif isinstance(stmt, Event):
                events.append(stmt)
        return events


class EidosDocument(object):
    def __init__(self, json_dict):
        self.tree = objectpath.Tree(json_dict)
        self.extractions = []
        self.sentences = {}
        self.entities = {}
        self.documents = {}
        self.coreferences = {}
        self.dct = None
        self._preprocess_extractions()

    def _preprocess_extractions(self):
        extractions = \
            self.tree.execute("$.extractions[(@.@type is 'Extraction')]")
        if not extractions:
            return
        # Listify for multiple reuse
        self.extractions = list(extractions)

        # Build a dictionary of entities
        entities = [e for e in self.extractions if 'Concept' in
                    e.get('labels', [])]
        self.entities = {entity['@id']: entity for entity in entities}

        # Build a dictionary of sentences and document creation times (DCTs)
        documents = self.tree.execute("$.documents[(@.@type is 'Document')]")
        self.sentences = {}
        for document in documents:
            title = document.get('title')
            self.documents[document['@id']] = {'title': title}
            # We stash the DCT here as a TimeContext object
            sentences = document.get('sentences', [])
            for sent in sentences:
                self.sentences[sent['@id']] = sent

        # Build a dictionary of coreferences
        for extraction in self.extractions:
            if 'Coreference' in extraction['labels']:
                reference = find_arg(extraction, 'reference')
                anchor = find_arg(extraction, 'anchor')
                self.coreferences[reference] = anchor


def _sanitize(text):
    """Return sanitized Eidos text field for human readability."""
    d = {'-LRB-': '(', '-RRB-': ')'}
    return re.sub('|'.join(d.keys()), lambda m: d[m.group(0)], text)


def find_arg(event, arg_type):
    """Return ID of the first argument of a given type"""
    obj_ids = find_args(event, arg_type)
    if not obj_ids:
        return None
    else:
        return obj_ids[0]


def find_args(event, arg_type):
    """Return IDs of all arguments of a given type"""
    args = event.get('arguments', {})
    obj_tags = [arg for arg in args if arg['type'] == arg_type]
    if obj_tags:
        return [o['value']['@id'] for o in obj_tags]
    else:
        return []
