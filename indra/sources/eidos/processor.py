from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import logging
import datetime
import objectpath
from indra.statements import Influence, Association, Concept, Evidence, \
    WorldContext, TimeContext, RefContext


logger = logging.getLogger(__name__)


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
        self.extractions = []
        self.sentences = {}
        self.entities = {}
        self.documents = {}
        self.coreferences = {}
        self.timexes = {}
        self.geolocs = {}
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
            dct = document.get('dct')
            title = document.get('title')
            self.documents[document['@id']] = {'title': title}
            # We stash the DCT here as a TimeContext object
            if dct is not None:
                self.dct = self.time_context_from_dct(dct)
                self.timexes[dct['@id']] = self.dct
            sentences = document.get('sentences', [])
            for sent in sentences:
                self.sentences[sent['@id']] = sent
                timexes = sent.get('timexes')
                if timexes:
                    for timex in timexes:
                        tc = self.time_context_from_timex(timex)
                        self.timexes[timex['@id']] = tc
                geolocs = sent.get('geolocs')
                if geolocs:
                    for geoloc in geolocs:
                        rc = self.ref_context_from_geoloc(geoloc)
                        self.geolocs[geoloc['@id']] = rc

        # Build a dictionary of coreferences
        for extraction in self.extractions:
            if 'Coreference' in extraction['labels']:
                reference = self.find_arg(extraction, 'reference')
                anchor = self.find_arg(extraction, 'anchor')
                self.coreferences[reference] = anchor

    def get_causal_relations(self):
        """Extract causal relations as Statements."""
        # Get the events that are labeled as directed and causal
        events = [e for e in self.extractions if
                  'DirectedRelation' in e['labels'] and
                  'Causal' in e['labels']]
        for event in events:
            # For now, just take the first source and first destination.
            # Later, might deal with hypergraph representation.
            subj_id = self.find_arg(event, 'source')
            obj_id = self.find_arg(event, 'destination')
            if subj_id is None or obj_id is None:
                continue

            # Resolve coreferences by ID
            subj_id = self.coreferences.get(subj_id, subj_id)
            obj_id = self.coreferences.get(obj_id, obj_id)

            # Get the actual entities
            subj = self.entities[subj_id]
            obj = self.entities[obj_id]

            subj_delta = self.extract_entity_states(subj.get('states', []))
            obj_delta = self.extract_entity_states(obj.get('states', []))

            evidence = self.get_evidence(event)

            # It is currently the case that time constraints and locations for
            #  concepts are better stored as annotations and the Evidence
            # level, we therefore move them over there.
            subj_timex = subj_delta.pop('time_context', None)
            obj_timex = obj_delta.pop('time_context', None)
            subj_geo = subj_delta.pop('geo_context', None)
            obj_geo = obj_delta.pop('geo_context', None)
            if subj_timex or subj_geo:
                wc = WorldContext(time=subj_timex,
                                  geo_location=subj_geo).to_json()
                evidence.annotations['subj_context'] = wc
            if obj_timex or obj_geo:
                wc = WorldContext(time=obj_timex,
                                  geo_location=obj_geo).to_json()
                evidence.annotations['obj_context'] = wc

            # In addition, for the time being we also put the adjectives and
            # polarities into annotations since they could otherwise get
            # squashed upon preassembly
            evidence.annotations['subj_adjectives'] = subj_delta['adjectives']
            evidence.annotations['obj_adjectives'] = obj_delta['adjectives']
            evidence.annotations['subj_polarity'] = subj_delta['polarity']
            evidence.annotations['obj_polarity'] = obj_delta['polarity']

            st = Influence(self.get_concept(subj), self.get_concept(obj),
                           subj_delta, obj_delta, evidence=[evidence])

            self.statements.append(st)

    def get_correlations(self):
        events = [e for e in self.extractions if
                  'UndirectedRelation' in e['labels'] and
                  'Correlation' in e['labels']]
        for event in events:
            # For now, just take the first source and first destination.
            # Later, might deal with hypergraph representation.
            arg_ids = self.find_args(event, 'argument')
            if len(arg_ids) != 2:
                logger.warning('Skipping correlation with not 2 arguments.')

            # Resolve coreferences by ID
            arg_ids = [self.coreferences.get(arg_id, arg_id)
                       for arg_id in arg_ids]

            # Get the actual entities
            args = [self.entities[arg_id] for arg_id in arg_ids]
            # Make Concepts from the entities
            members = [self.get_concept(arg) for arg in args]
            # Get the evidence
            evidence = self.get_evidence(event)

            st = Association(members, evidence=[evidence])
            self.statements.append(st)

    def get_evidence(self, event):
        """Return the Evidence object for the INDRA Statment."""
        provenance = event.get('provenance')

        # First try looking up the full sentence through provenance
        text = None
        context = None
        if provenance:
            sentence_tag = provenance[0].get('sentence')
            if sentence_tag and '@id' in sentence_tag:
                sentence_id = sentence_tag['@id']
                sentence = self.sentences.get(sentence_id)
                if sentence is not None:
                    text = _sanitize(sentence['text'])
                # Get temporal constraints if available
                timexes = sentence.get('timexes', [])
                if timexes:
                    # We currently handle just one timex per statement
                    timex = timexes[0]
                    tc = self.time_context_from_timex(timex)
                    context = WorldContext(time=tc)
                # Get geolocation if available
                geolocs = sentence.get('geolocs', [])
                if geolocs:
                    geoloc = geolocs[0]
                    rc = self.ref_context_from_geoloc(geoloc)
                    if context:
                        context.geo_location = rc
                    else:
                        context = WorldContext(geo_location=rc)

            # Here we try to get the title of the document and set it
            # in the provenance
            doc_id = provenance[0].get('document', {}).get('@id')
            if doc_id:
                title = self.documents.get(doc_id, {}).get('title')
                if title:
                    provenance[0]['document']['title'] = title

        annotations = {'found_by': event.get('rule'),
                       'provenance': provenance}
        if self.dct is not None:
            annotations['document_creation_time'] = self.dct.to_json()

        epistemics = {}
        negations = self.get_negation(event)
        hedgings = self.get_hedging(event)
        if hedgings:
            epistemics['hedgings'] = hedgings
        if negations:
            # This is the INDRA standard to show negation
            epistemics['negated'] = True
            # But we can also save the texts associated with the negation
            # under annotations, just in case it's needed
            annotations['negated_texts'] = negations

        # If that fails, we can still get the text of the event
        if text is None:
            text = _sanitize(event.get('text'))

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
        time_context = None
        geo_context = None
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
            if state['type'] == 'TIMEX':
                time_context = self.time_context_from_ref(state)
            elif state['type'] == 'LocationExp':
                geo_context = self.geo_context_from_ref(state)
        return {'polarity': polarity, 'adjectives': adjectives,
                'time_context': time_context, 'geo_context': geo_context}

    @staticmethod
    def get_groundings(entity):
        """Return groundings as db_refs for an entity."""
        def get_grounding_entries(grounding):
            if not grounding:
                return None

            entries = []
            values = grounding.get('values', [])
            # Values could still have been a None entry here
            if values:
                for entry in values:
                    ont_concept = entry.get('ontologyConcept')
                    value = entry.get('value')
                    if ont_concept is None or value is None:
                        continue
                    entries.append((ont_concept, value))
            return entries

        # Save raw text and Eidos scored groundings as db_refs
        db_refs = {'TEXT': entity['text']}
        groundings = entity.get('groundings')
        if not groundings:
            return db_refs
        for g in groundings:
            entries = get_grounding_entries(g)
            # Only add these groundings if there are actual values listed
            if entries:
                key = g['name'].upper()
                if key == 'UN':
                    db_refs[key] = [(s[0].replace(' ', '_'), s[1])
                                    for s in entries]
                else:
                    db_refs[key] = entries
        return db_refs

    @staticmethod
    def get_concept(entity):
        """Return Concept from an Eidos entity."""
        # Use the canonical name as the name of the Concept
        name = entity['canonicalName']
        db_refs = EidosProcessor.get_groundings(entity)
        concept = Concept(name, db_refs=db_refs)
        return concept

    @staticmethod
    def find_arg(event, arg_type):
        """Return ID of the first argument of a given type"""
        obj_ids = EidosProcessor.find_args(event, arg_type)
        if not obj_ids:
            return None
        else:
            return obj_ids[0]

    @staticmethod
    def find_args(event, arg_type):
        """Return IDs of all arguments of a given type"""
        args = event.get('arguments', {})
        obj_tags = [arg for arg in args if arg['type'] == arg_type]
        if obj_tags:
            return [o['value']['@id'] for o in obj_tags]
        else:
            return []

    @staticmethod
    def time_context_from_dct(dct):
        """Return a time context object given a DCT entry."""
        time_text = dct.get('text')
        start = _get_time_stamp(dct.get('start'))
        end = _get_time_stamp(dct.get('end'))
        duration = dct.get('duration')
        tc = TimeContext(text=time_text, start=start, end=end,
                         duration=duration)
        return tc

    def time_context_from_ref(self, timex):
        """Return a time context object given a timex reference entry."""
        # If the timex has a value set, it means that it refers to a DCT or
        # a TimeExpression e.g. "value": {"@id": "_:DCT_1"} and the parameters
        # need to be taken from there
        value = timex.get('value')
        if value:
            # Here we get the TimeContext directly from the stashed DCT
            # dictionary
            tc = self.timexes.get(value['@id'])
            return tc
        return None

    def geo_context_from_ref(self, ref):
        """Return a ref context object given a location reference entry."""
        value = ref.get('value')
        if value:
            # Here we get the RefContext from the stashed geoloc dictionary
            rc = self.geolocs.get(value['@id'])
            return rc
        return None

    @staticmethod
    def time_context_from_timex(timex):
        """Return a TimeContext object given a timex entry."""
        time_text = timex.get('text')
        constraint = timex['intervals'][0]
        start = _get_time_stamp(constraint.get('start'))
        end = _get_time_stamp(constraint.get('end'))
        duration = constraint['duration']
        tc = TimeContext(text=time_text, start=start, end=end,
                         duration=duration)
        return tc

    @staticmethod
    def ref_context_from_geoloc(geoloc):
        """Return a RefContext object given a geoloc entry."""
        text = geoloc.get('text')
        geoid = geoloc.get('geoID')
        rc = RefContext(name=text, db_refs={'GEOID': geoid})
        return rc


def _sanitize(text):
    """Return sanitized Eidos text field for human readability."""
    d = {'-LRB-': '(', '-RRB-': ')'}
    return re.sub('|'.join(d.keys()), lambda m: d[m.group(0)], text)


def _get_time_stamp(entry):
    """Return datetime object from a timex constraint start/end entry.

    Example string format to convert: 2018-01-01T00:00
    """
    if not entry or entry == 'Undef':
        return None
    try:
        dt = datetime.datetime.strptime(entry, '%Y-%m-%dT%H:%M')
    except Exception as e:
        logger.debug('Could not parse %s format' % entry)
        return None
    return dt
