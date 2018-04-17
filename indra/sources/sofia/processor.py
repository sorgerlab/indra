from indra.statements import Influence, Concept, Evidence

pos_rels = ['provide', 'led', 'lead', 'driv', 'support', 'enabl', 'develop']
neg_rels = ['restrict', 'worsen', 'declin', 'limit', 'constrain',
            'decreas', 'hinder', 'deplet', 'reduce', 'hamper']
neu_rels = ['affect']


class SofiaProcessor(object):
    def __init__(self, rows):
        self.statements = []
        header = [cell.value for cell in next(rows)]
        for row in rows:
            stmt = self._process_row(header, [r.value for r in row])
            if stmt:
                self.statements.append(stmt)

    @staticmethod
    def _process_row(header, row):
        row_dict = {h: v for h, v in zip(header, row)}
        subj = row_dict.get('Agent')
        obj = row_dict.get('Patient')
        if not obj or not subj:
            return None
        rel = row_dict.get('Relation')
        if _in_rels(rel, pos_rels):
            pol = 1
        elif _in_rels(rel, neg_rels):
            pol = -1
        elif _in_rels(rel, neu_rels):
            pol = None
        else:
            return None

        subj_concept = Concept(subj)
        obj_concept = Concept(obj)
        text = row_dict.get('Sentence')
        annot_keys = ['Relation', 'Event_Type', 'Location', 'Time']
        annots = {k: row_dict.get(k) for k in annot_keys}
        ev = Evidence(source_api='sofia', annotations=annots, text=text)
        stmt = Influence(subj_concept, obj_concept, evidence=[ev])
        stmt.obj_delta['polarity'] = pol
        return stmt


def _in_rels(value, rels):
    for rel in rels:
        if value.startswith(rel):
            return True
    return False
