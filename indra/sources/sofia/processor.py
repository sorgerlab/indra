import itertools
from indra.statements import Influence, Concept, Evidence, WorldContext, \
    TimeContext, RefContext

pos_rels = ['provide', 'led', 'lead', 'driv', 'support', 'enabl', 'develop']
neg_rels = ['restrict', 'worsen', 'declin', 'limit', 'constrain',
            'decreas', 'hinder', 'deplet', 'reduce', 'hamper']
neu_rels = ['affect', 'impact', 'due', 'caus', 'because']


class SofiaProcessor(object):
    def __init__(self, input_object):
        self.statements = []
        self._events = {}
        if isinstance(input_object, list):  # Sofia API
            self._events = self._process_events_json(input_object)
            self.statements = self._process_relations_json(input_object)
        elif isinstance(input_object, tuple):  # excel
            relation_rows, event_rows, entity_rows = input_object
            self._events = self._process_events_excel(event_rows)
            self.statements = self._process_relations_excel(relation_rows)

    def _process_events_json(self, json_list):
        event_dict = {}
        for _dict in json_list:
            events = _dict['Events']
            for event in events:
                # subj_entries = event.get('Agent')
                # object_entries = event.get('Patient')
                event_index = event.get('Event Index')
                event_dict[event_index] = self._process_event(event)

        return event_dict

    def _process_relations_json(self, json_list):
        stmts = []
        for _dict in json_list:
            json_relations = _dict['Causal']
            for rel_dict in json_relations:
                stmt_list = self._build_stmt(rel_dict)
                if not stmt_list:
                    continue
                stmts = stmts + stmt_list
        return stmts

    def _process_events_excel(self, event_rows):
        header = [cell.value for cell in next(event_rows)]
        event_dict = {}
        for row in event_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(header, row_values)}
            # subj_entries = row_dict.get('Agent')
            # obj_entries = row_dict.get('Patient')
            event_index = row_dict.get('Event Index')
            event_dict[event_index] = self._process_event(row_dict)
        return event_dict

    def _process_relations_excel(self, relation_rows):
        header = [cell.value for cell in next(relation_rows)]
        stmts = []
        for row in relation_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(header, row_values)}
            stmt_list = self._build_stmt(row_dict)
            if not stmt_list:
                continue
            stmts = stmts + stmt_list
        return stmts

    @staticmethod
    def _process_event(_dict):
        return {'Event_Type': _dict.get('Event_Type'),
                'Relation': _dict.get('Relation'),
                'Location': _dict.get('Location'),
                'Time': _dict.get('Time')}

    def _build_stmt(self, rel_dict):
        stmt_list = []
        cause_entries = rel_dict.get('Cause Index')
        effect_entries = rel_dict.get('Effect Index')

        # FIXME: Handle cases in which there is a missing cause/effect
        if not cause_entries or not effect_entries:
            return []
        causes = [c.strip() for c in cause_entries.split(',')]
        effects = [e.strip() for e in effect_entries.split(',')]
        rel = rel_dict.get('Relation')
        if _in_rels(rel, pos_rels):
            pol = 1
        elif _in_rels(rel, neg_rels):
            pol = -1
        elif _in_rels(rel, neu_rels):
            pol = None
        # If we don't recognize this relation, we don't get any
        # statements
        else:
            return []

        text = rel_dict.get('Sentence')
        annot_keys = ['Relation']
        annots = {k: rel_dict.get(k) for k in annot_keys}
        ref = rel_dict.get('Source_File')

        for cause_idx, effect_idx in itertools.product(causes,
                                                       effects):
            cause_name = self._events[cause_idx]['Relation']
            cause_grounding = self._events[cause_idx]['Event_Type']
            effect_name = self._events[effect_idx]['Relation']
            effect_grounding = self._events[effect_idx]['Event_Type']
            cause_concept = Concept(cause_name,
                                    db_refs={'TEXT': cause_name})
            if cause_grounding:
                cause_concept.db_refs['SOFIA'] = cause_grounding
            effect_concept = Concept(effect_name,
                                     db_refs={'TEXT': effect_name})
            if effect_grounding:
                effect_concept.db_refs['SOFIA'] = effect_grounding

            # NOTE: Extract context. The basic issue is that
            # time/location
            # here is given at the event level, not at the relation
            # level, and so we need to choose which event's context
            # we will associate with the relation
            def choose_context(context_type):
                locs = [self._events[cause_idx].get(context_type),
                        self._events[effect_idx].get(context_type)]
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

            ev = Evidence(source_api='sofia', pmid=ref,
                          annotations=annots, text=text,
                          context=context)
            stmt = Influence(cause_concept, effect_concept,
                             evidence=[ev])
            # Assume unknown polarity on the subject, put the overall
            # polarity in the sign of the object
            stmt.subj_delta['polarity'] = None
            stmt.obj_delta['polarity'] = pol

            stmt_list.append(stmt)
        return stmt_list


def _in_rels(value, rels):
    for rel in rels:
        if value.lower().startswith(rel):
            return True
    return False
