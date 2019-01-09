import itertools

from indra.statements import Influence, Concept, Evidence, WorldContext, \
    TimeContext, RefContext

pos_rels = ['provide', 'led', 'lead', 'driv', 'support', 'enabl', 'develop']
neg_rels = ['restrict', 'worsen', 'declin', 'limit', 'constrain',
            'decreas', 'hinder', 'deplet', 'reduce', 'hamper']
neu_rels = ['affect', 'impact', 'due', 'caus', 'because']


class SofiaJsonProcessor(object):
    def __init__(self, json_events, json_relations, json_entities):
        self.statements = []
        self._events = self._process_events(json_events)
        self.statements = self._process_relations(json_relations, self._events)

    @staticmethod
    def _process_events(json_events):
        event_dict = {}
        for event in json_events:
            subj_entries = event.get('Agent')
            object_entries = event.get('Patient')
            event_index = event.get('Event Index')
            event_type = event.get('Event_Type')
            relation = event.get('Relation')
            location = event.get('Location')
            time = event.get('Time')
            event_dict[event_index] = {'Event_Type': event_type,
                                       'Relation': relation,
                                       'Location': location,
                                       'Time': time}
        return event_dict

    @staticmethod
    def _process_relations(json_relations, event_dict):
        stmts = []
        for rel_dict in json_relations:
            cause_entries = rel_dict.get('Cause Index')
            effect_entries = rel_dict.get('Effect Index')
            if not cause_entries or not effect_entries:
                continue
            causes = [c.strip() for c in cause_entries.split(',')]
            effects = [e.strip() for e in effect_entries.split(',')]
            rel = rel_dict.get('Relation')
            if _in_rels(rel, pos_rels):
                pol = 1
            elif _in_rels(rel, neg_rels):
                pol = -1
            elif _in_rels(rel, neu_rels):
                pol = None
            # If we don't recognize this relation, we don't get any statements
            else:
                continue

            text = rel_dict.get('Sentence')
            annot_keys = ['Relation']
            annots = {k: rel_dict.get(k) for k in annot_keys}
            if __name__ == '__main__':
                ref = rel_dict.get('Source_File')

            all_contexts = []
            for cause_idx, effect_idx in itertools.product(causes, effects):
                cause_name = event_dict[cause_idx]['Relation']
                cause_grounding = event_dict[cause_idx]['Event_Type']
                effect_name = event_dict[effect_idx]['Relation']
                effect_grounding = event_dict[effect_idx]['Event_Type']
                cause_concept = Concept(cause_name,
                                        db_refs={'TEXT': cause_name})
                if cause_grounding:
                    cause_concept.db_refs['SOFIA'] = cause_grounding
                effect_concept = Concept(effect_name,
                                         db_refs={'TEXT': effect_name})
                if effect_grounding:
                    effect_concept.db_refs['SOFIA'] = effect_grounding

                # NOTE: Extract context. The basic issue is that time/location
                # here is given at the event level, not at the relation
                # level, and so we need to choose which event's context
                # we will associate with the relation
                def choose_context(context_type):
                    locs = [event_dict[cause_idx].get(context_type),
                            event_dict[effect_idx].get(context_type)]
                    if locs[0]:
                        return locs[0].strip()
                    elif locs[1]:
                        return locs[1].strip()
                    else:
                        return None
                context = WorldContext()
                location = choose_context('Location')
                if location:
                    context.location = RefContext(name=location)
                time = choose_context('Time')
                if time:
                    context.time = TimeContext(text=time)
                # Overwrite blank context
                if not context:
                    context = None

                ev = Evidence(source_api='sofia', pmid=ref, annotations=annots,
                              text=text, context=context)
                stmt = Influence(cause_concept, effect_concept, evidence=[ev])
                # Assume unknown polarity on the subject, put the overall
                # polarity in the sign of the object
                stmt.subj_delta['polarity'] = None
                stmt.obj_delta['polarity'] = pol

                stmts.append(stmt)
        return stmts


class SofiaProcessor(object):
    def __init__(self, relation_rows, event_rows, entity_rows):
        self.statements = []
        self._events = self._process_events(event_rows)
        self.statements = self._process_relations(relation_rows, self._events)

    @staticmethod
    def _process_events(event_rows):
        header = [cell.value for cell in next(event_rows)]
        event_dict = {}
        for row in event_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(header, row_values)}
            subj_entries = row_dict.get('Agent')
            obj_entries = row_dict.get('Patient')
            relation = row_dict.get('Relation')
            event_type = row_dict.get('Event_Type')
            event_index = row_dict.get('Event Index')
            location = row_dict.get('Location')
            time = row_dict.get('Time')
            event_dict[event_index] = {'Event_Type': event_type,
                                       'Relation': relation,
                                       'Location': location,
                                       'Time': time}
        return event_dict

    @staticmethod
    def _process_relations(relation_rows, event_dict):
        header = [cell.value for cell in next(relation_rows)]
        stmts = []
        for row in relation_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(header, row_values)}
            cause_entries = row_dict.get('Cause Index')
            effect_entries = row_dict.get('Effect Index')

            # FIXME: Handle cases in which there is a missing cause/effect
            if not cause_entries or not effect_entries:
                continue
            causes = [c.strip() for c in cause_entries.split(',')]
            effects = [e.strip() for e in effect_entries.split(',')]

            rel = row_dict.get('Relation')
            if _in_rels(rel, pos_rels):
                pol = 1
            elif _in_rels(rel, neg_rels):
                pol = -1
            elif _in_rels(rel, neu_rels):
                pol = None
            # If we don't recognize this relation, we don't get any statements
            else:
                continue

            text = row_dict.get('Sentence')
            #annot_keys = ['Relation', 'Event_Type', 'Location', 'Time']
            #annots = {k: row_dict.get(k) for k in annot_keys}
            annot_keys = ['Relation']
            annots = {k: row_dict.get(k) for k in annot_keys}
            ref = row_dict.get('Source_File')

            all_contexts = []
            for cause_idx, effect_idx in itertools.product(causes, effects):
                cause_name = event_dict[cause_idx]['Relation']
                cause_grounding = event_dict[cause_idx]['Event_Type']
                effect_name = event_dict[effect_idx]['Relation']
                effect_grounding = event_dict[effect_idx]['Event_Type']
                cause_concept = Concept(cause_name,
                                        db_refs={'TEXT': cause_name})
                if cause_grounding:
                    cause_concept.db_refs['SOFIA'] = cause_grounding
                effect_concept = Concept(effect_name,
                                         db_refs={'TEXT': effect_name})
                if effect_grounding:
                    effect_concept.db_refs['SOFIA'] = effect_grounding

                # NOTE: Extract context. The basic issue is that time/location
                # here is given at the event level, not at the relation
                # level, and so we need to choose which event's context
                # we will associate with the relation
                def choose_context(context_type):
                    locs = [event_dict[cause_idx].get(context_type),
                            event_dict[effect_idx].get(context_type)]
                    if locs[0]:
                        return locs[0].strip()
                    elif locs[1]:
                        return locs[1].strip()
                    else:
                        return None
                context = WorldContext()
                location = choose_context('Location')
                if location:
                    context.location = RefContext(name=location)
                time = choose_context('Time')
                if time:
                    context.time = TimeContext(text=time)
                # Overwrite blank context
                if not context:
                    context = None

                ev = Evidence(source_api='sofia', pmid=ref, annotations=annots,
                              text=text, context=context)
                stmt = Influence(cause_concept, effect_concept, evidence=[ev])
                # Assume unknown polarity on the subject, put the overall
                # polarity in the sign of the object
                stmt.subj_delta['polarity'] = None
                stmt.obj_delta['polarity'] = pol

                stmts.append(stmt)
        return stmts


def _in_rels(value, rels):
    for rel in rels:
        if value.lower().startswith(rel):
            return True
    return False
