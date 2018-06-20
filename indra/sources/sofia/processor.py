import itertools
from indra.statements import Influence, Concept, Evidence

pos_rels = ['provide', 'led', 'lead', 'driv', 'support', 'enabl', 'develop']
neg_rels = ['restrict', 'worsen', 'declin', 'limit', 'constrain',
            'decreas', 'hinder', 'deplet', 'reduce', 'hamper']
neu_rels = ['affect', 'impact', 'due', 'caus', 'because']


class SofiaProcessor(object):
    def __init__(self, rows):
        self.statements = []
        header = [cell.value for cell in next(rows)]
        for row in rows:
            stmts = self._process_row(header, [r.value for r in row])
            if stmts:
                self.statements.extend(stmts)

    @staticmethod
    def _process_row(header, row):
        row_dict = {h: v for h, v in zip(header, row)}
        cause_entries = row_dict.get('Cause')
        effect_entries = row_dict.get('Effect')
        # FIXME: Handle cases in which there is a missing cause/effect
        if not cause_entries or not effect_entries:
            return None
        causes = [c.strip() for c in cause_entries.split(',')]
        effects = [e.strip() for e in effect_entries.split(',')]

        #subj = row_dict.get('Agent')
        #obj = row_dict.get('Patient')
        rel = row_dict.get('Relation')
        if _in_rels(rel, pos_rels):
            pol = 1
        elif _in_rels(rel, neg_rels):
            pol = -1
        elif _in_rels(rel, neu_rels):
            pol = None
        # If we don't recognize this relation, we don't get any statements
        else:
            return []

        text = row_dict.get('Sentence')
        #annot_keys = ['Relation', 'Event_Type', 'Location', 'Time']
        #annots = {k: row_dict.get(k) for k in annot_keys}
        annot_keys = ['Relation']
        annots = {k: row_dict.get(k) for k in annot_keys}
        ref = row_dict.get('Source_File')
        ev = Evidence(source_api='sofia', pmid=ref, annotations=annots,
                      text=text)

        stmts = []
        for cause, effect in itertools.product(causes, effects):
            cause_concept = Concept(cause, db_refs={'TEXT': cause,
                                                    'SOFIA': cause})
            effect_concept = Concept(effect, db_refs={'TEXT': effect,
                                                      'SOFIA': effect})
            stmt = Influence(cause_concept, effect_concept, evidence=[ev])
            stmt.obj_delta['polarity'] = pol
            stmts.append(stmt)

        return stmts



def _in_rels(value, rels):
    for rel in rels:
        if value.lower().startswith(rel):
            return True
    return False
