import itertools
from indra.statements import Influence, Concept, Evidence

pos_rels = ['provide', 'led', 'lead', 'driv', 'support', 'enabl', 'develop']
neg_rels = ['restrict', 'worsen', 'declin', 'limit', 'constrain',
            'decreas', 'hinder', 'deplet', 'reduce', 'hamper']
neu_rels = ['affect', 'impact', 'due', 'caus', 'because']


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
            event_dict[event_index] = {'Event_Type': event_type,
                                       'Relation': relation}
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
            ev = Evidence(source_api='sofia', pmid=ref, annotations=annots,
                          text=text)

            for cause_index, effect_index in itertools.product(causes, effects):
                cause_name = event_dict[cause_index]['Relation']
                cause_grounding = event_dict[cause_index]['Event_Type']
                effect_name = event_dict[effect_index]['Relation']
                effect_grounding = event_dict[effect_index]['Event_Type']
                cause_concept = Concept(cause_name,
                                        db_refs={'TEXT': cause_name,
                                                 'SOFIA': cause_grounding})
                effect_concept = Concept(effect_name,
                                         db_refs={'TEXT': effect_name,
                                                  'SOFIA': effect_grounding})
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
