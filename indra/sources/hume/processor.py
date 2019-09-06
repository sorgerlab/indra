import os
import rdflib
import logging
import objectpath
from datetime import datetime, timedelta


from indra.statements import Concept, Event, Influence, TimeContext, \
    RefContext, WorldContext, Evidence, QualitativeDelta, MovementContext, \
    Migration, QuantitativeState


logger = logging.getLogger(__name__)


# List out relation types and their default (implied) polarities.
polarities = {'causation': 1, 'precondition': 1, 'catalyst': 1,
              'mitigation': -1, 'prevention': -1,
              'temporallyPrecedes': None}


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
        self._get_documents()
        self.relation_subj_obj_ids = []

    def extract_relations(self):
        relations = self._find_relations()
        for relation_type, relation in relations:
            # Extract concepts and contexts.
            subj = self._get_event_and_context(relation, arg_type='source')
            obj = self._get_event_and_context(relation, arg_type='destination')

            if not subj.concept or not obj.concept:
                continue

            # Apply the naive polarity from the type of statement. For the
            # purpose of the multiplication here, if obj.delta.polarity is
            # None to begin with, we assume it is positive
            obj_pol = obj.delta.polarity
            obj_pol = obj_pol if obj_pol is not None else 1
            rel_pol = polarities[relation_type]
            obj.delta.polarity = rel_pol * obj_pol if rel_pol else None

            evidence = self._get_evidence(relation, get_states(relation))

            st = Influence(subj, obj, evidence=evidence)
            self.eid_stmt_dict[relation['@id']] = st
            self.statements.append(st)

    def extract_events(self):
        events = self._find_events()
        for event in events:
            evidence = self._get_evidence(event, get_states(event))
            stmt = self._get_event_and_context(event, eid=event['@id'],
                                               evidence=evidence)
            self.eid_stmt_dict[event['@id']] = stmt
            self.statements.append(stmt)

    def _find_events(self):
        """Find standalone events and return them in a list."""
        # First populate self.concept_dict and self.relations_subj_obj_ids
        if not self.relation_dict or not self.concept_dict or \
                not self.relation_subj_obj_ids:
            self._find_relations()

        # Check if events are part of relations
        events = []
        for e in self.concept_dict.values():
            label_set = set(e.get('labels', []))
            if 'Event' in label_set:
                if e['@id'] not in self.relation_subj_obj_ids:
                    events.append(e)

        if not events:
            logger.info('No standalone events found.')
        else:
            logger.info('%d standalone events found.' % len(events))

        return events

    def _find_relations(self):
        """Find all relevant relation elements and return them in a list."""
        # Get all extractions
        extractions = \
            list(self.tree.execute("$.extractions[(@.@type is 'Extraction')]"))

        # Get relations from extractions
        relations = []
        for e in extractions:
            label_set = set(e.get('labels', []))
            # If this is a DirectedRelation
            if 'DirectedRelation' in label_set:
                self.relation_dict[e['@id']] = e
                subtype = e.get('subtype')
                if any(t in subtype for t in polarities.keys()):
                    relations.append((subtype, e))
                    # Save IDs of relation's subject and object
                    if e['arguments']:
                        for a in e['arguments']:
                            if a['type'] == 'source' or \
                                    a['type'] == 'destination':
                                self.relation_subj_obj_ids.append(
                                    a['value']['@id'])
            # If this is an Event or an Entity
            if {'Event', 'Entity'} & label_set:
                self.concept_dict[e['@id']] = e

        if not relations and not self.relation_dict:
            logger.info("No relations found.")
        else:
            logger.info('%d relations of types %s found'
                        % (len(relations), ', '.join(polarities.keys())))
            logger.info('%d relations in dict.' % len(self.relation_dict))
            logger.info('%d concepts found.' % len(self.concept_dict))
        return relations

    def _get_documents(self):
        """Populate sentences attribute with a dict keyed by document id."""
        documents = self.tree.execute("$.documents")
        for doc in documents:
            sentences = {s['@id']: s['text'] for s in doc.get('sentences', [])}
            self.document_dict[doc['@id']] = {'sentences': sentences,
                                              'location': doc['location']}

    def _make_world_context(self, entity):
        """Get place and time info from the json for this entity."""
        loc_context = None
        time_context = None

        # Look for time and place contexts.
        for argument in entity["arguments"]:
            if argument["type"] in {"has_location", "has_origin_location",
                                    "has_destination_location",
                                    "has_intermediate_location"}:
                entity_id = argument["value"]["@id"]
                loc_entity = self.concept_dict[entity_id]
                loc_context = _resolve_geo(loc_entity)
            if argument["type"] in {"has_time", "has_start_time",
                                    "has_end_time"}:
                entity_id = argument["value"]["@id"]
                temporal_entity = self.concept_dict[entity_id]
                time_context = _resolve_time(temporal_entity)

        # Put context together
        context = None
        if loc_context or time_context:
            context = WorldContext(time=time_context, geo_location=loc_context)

        return context

    def _make_movement_context(self, entity):
        movement_locations = list()
        time_context = None
        # Use None for quantitative_state if no information found, default
        # value will be assigned when creating a Statement
        quantitative_state = None
        for argument in entity['arguments']:
            entity_id = argument["value"]["@id"]
            hume_entity = self.concept_dict[entity_id]
            if argument['type'] in {"has_actor", "has_affected_actor",
                                    "has_active_actor"}:
                for count in hume_entity.get('counts', list()):
                    quantitative_state = QuantitativeState(
                        entity="person", value=count['value'],
                        unit=count['unit'], modifier=count['modifier'])
            if argument['type'] == "has_origin_location":
                movement_locations.append(
                    {'location': _resolve_geo(hume_entity), 'role': 'origin'})
            if argument['type'] == 'has_destination_location':
                movement_locations.append(
                    {'location': _resolve_geo(hume_entity),
                     'role': 'destination'})
            if argument['type'] in {"has_time", "has_start_time",
                                    "has_end_time"}:
                time_context = _resolve_time(hume_entity)
        return MovementContext(locations=movement_locations,
                               time=time_context), quantitative_state

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
        db_refs = _get_grounding(entity)
        concept = Concept(name, db_refs=db_refs)
        metadata = {arg['type']: arg['value']['@id']
                    for arg in entity['arguments']}

        return concept, metadata

    def _get_event_and_context(self, event, eid=None, arg_type=None,
                               evidence=None):
        """Return an INDRA Event based on an event entry."""
        if not eid:
            eid = _choose_id(event, arg_type)
        ev = self.concept_dict[eid]
        concept, metadata = self._make_concept(ev)

        is_migration_event = False
        hume_grounding = {x[0] for x in concept.db_refs['WM']}
        for grounding_en in hume_grounding:
            if "wm/concept/causal_factor/social_and_political/migration" in \
                    grounding_en:
                is_migration_event = True
        if is_migration_event:
            movement_context, quantitative_state = (
                self._make_movement_context(ev))
            event_obj = Migration(concept, delta=quantitative_state,
                                  context=movement_context, evidence=evidence)
        else:
            ev_delta = QualitativeDelta(
                polarity=get_polarity(ev), adjectives=None)
            context = self._make_world_context(ev)
            event_obj = Event(concept, delta=ev_delta, context=context,
                              evidence=evidence)
        return event_obj

    def _get_evidence(self, event, adjectives):
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
            'event_type': os.path.basename(event.get('type')),
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


def _get_grounding(entity):
    """Return Hume grounding."""
    db_refs = {'TEXT': entity['text']}
    groundings = entity.get('grounding')
    if not groundings:
        return db_refs
    # Get rid of leading slash
    groundings = [(x['ontologyConcept'][1:], x['value']) for x in groundings]
    grounding_entries = sorted(list(set(groundings)),
                               key=lambda x: (x[1], x[0].count('/'), x[0]),
                               reverse=True)
    # We could get an empty list here in which case we don't add the
    # grounding
    if grounding_entries:
        db_refs['WM'] = grounding_entries
    return db_refs


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


def _resolve_geo(hume_loc_entity):
    place = hume_loc_entity.get('canonicalName', hume_loc_entity.get('text'))
    geo_id = hume_loc_entity.get('geoname_id', None)
    if geo_id is not None:
        return RefContext(name=place, db_refs={"GEOID": geo_id})
    else:
        return RefContext(place)


def _resolve_time(hume_temporal_entity):
    text = hume_temporal_entity['mentions'][0]['text']
    if len(hume_temporal_entity.get("timeInterval", [])) < 1:
        return TimeContext(text=text)
    time = hume_temporal_entity["timeInterval"][0]
    start = datetime.strptime(time['start'], '%Y-%m-%dT%H:%M')
    end = datetime.strptime(time['end'], '%Y-%m-%dT%H:%M')
    end = end + timedelta(minutes=1)
    duration = int(time['duration'])
    return TimeContext(text=text, start=start, end=end,
                       duration=duration)
