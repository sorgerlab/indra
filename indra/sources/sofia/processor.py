import itertools
from indra.statements import Influence, Concept, Event, Evidence, \
    WorldContext, TimeContext, RefContext

pos_rels = ['provide', 'led', 'lead', 'driv', 'support', 'enabl', 'develop']
neg_rels = ['restrict', 'worsen', 'declin', 'limit', 'constrain',
            'decreas', 'hinder', 'deplet', 'reduce', 'hamper']
neu_rels = ['affect', 'impact', 'due', 'caus', 'because']


class SofiaProcessor(object):
    @staticmethod
    def process_event(event_dict):
        return {'Event_Type': event_dict.get('Event_Type'),
                'Relation': event_dict.get('Relation'),
                'Location': event_dict.get('Location'),
                'Time': event_dict.get('Time')}

    def _build_influences(self, rel_dict):
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

        for cause_idx, effect_idx in itertools.product(causes, effects):
            subj = self.get_event(self._events[cause_idx])
            obj = self.get_event(self._events[effect_idx])

            ev = Evidence(source_api='sofia', pmid=ref,
                          annotations=annots, text=text)
            stmt = Influence(subj, obj, evidence=[ev])
            # Assume unknown polarity on the subject, put the overall
            # polarity in the sign of the object
            stmt.subj.delta['polarity'] = None
            stmt.obj.delta['polarity'] = pol

            stmt_list.append(stmt)
        return stmt_list

    @staticmethod
    def get_event(event_entry):
        name = event_entry['Relation']
        concept = Concept(name, db_refs={'TEXT': name})
        grounding = event_entry['Event_Type']
        if grounding:
            concept.db_refs['SOFIA'] = grounding
        context = WorldContext()
        time = event_entry.get('Time')
        if time:
            context.time = TimeContext(text=time.strip())
        loc = event_entry.get('Location')
        if loc:
            context.geo_location = RefContext(name=loc)
        event = Event(concept, context=context)
        return event


class SofiaJsonProcessor(SofiaProcessor):
    def __init__(self, json_list):
        self._events = self.process_events(json_list)
        self.statements = self.process_relations(json_list)

    def process_events(self, json_list):
        event_dict = {}
        for _dict in json_list:
            events = _dict['Events']
            for event in events:
                event_index = event.get('Event Index')
                event_dict[event_index] = self.process_event(event)

        return event_dict

    def process_relations(self, json_list):
        stmts = []
        for _dict in json_list:
            json_relations = _dict['Causal']
            for rel_dict in json_relations:
                stmt_list = self._build_influences(rel_dict)
                if not stmt_list:
                    continue
                stmts = stmts + stmt_list
        return stmts


class SofiaExcelProcessor(SofiaProcessor):
    def __init__(self, event_rows):
        self._events = self.process_events(event_rows)
        self.statements = []

    def process_events(self, event_rows):
        header = [cell.value for cell in next(event_rows)]
        event_dict = {}
        for row in event_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(header, row_values)}
            event_index = row_dict.get('Event Index')
            event_dict[event_index] = self.process_event(row_dict)
        return event_dict

    def extract_relations(self, relation_rows):
        header = [cell.value for cell in next(relation_rows)]
        stmts = []
        for row in relation_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(header, row_values)}
            stmt_list = self._build_influences(row_dict)
            if not stmt_list:
                continue
            stmts = stmts + stmt_list
        for stmt in stmts:
            self.statements.append(stmt)

    def extract_events(self, event_rows, relation_rows):
        # First find indexes of event that are part of causal relations
        relation_events = []
        rel_header = [cell.value for cell in next(relation_rows)]
        for row in relation_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(rel_header, row_values)}
            cause_entries = row_dict.get('Cause Index')
            effect_entries = row_dict.get('Effect Index')
            causes = [c.strip() for c in cause_entries.split(',')]
            effects = [e.strip() for e in effect_entries.split(',')]
            for ci in causes:
                relation_events.append(ci)
            for ei in effects:
                relation_events.append(ei)
        event_header = [cell.value for cell in next(event_rows)]
        row_count = 0
        for row in event_rows:
            row_count += 1
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(event_header, row_values)}
            event_index = row_dict.get('Event Index')
            if event_index in relation_events:
                continue
            event = self.get_event(self._events[event_index])
            self.statements.append(event)


def _in_rels(value, rels):
    for rel in rels:
        if value.lower().startswith(rel):
            return True
    return False
