import openpyxl
from indra.statements import Influence, Concept, Evidence

pos_rels = ['provide', 'led', 'lead', 'driv', 'support', 'enabl', 'develop']
neg_rels = ['restrict', 'worsen', 'declin', 'limit', 'constrain',
            'decreas', 'hinder', 'deplet', 'reduce', 'hamper']
neu_rels = ['affect']

def process_table(fname, sheet_name):
    book = openpyxl.load_workbook(fname, read_only=True)
    sheet = book[sheet_name]
    rows = sheet.rows
    header = [cell.value for cell in next(rows)]

    stmts = []
    for row in rows:
        stmt = _process_row(header, [r.value for r in row])
        if stmt:
            stmts.append(stmt)
    return stmts


def _in_rels(value, rels):
    for rel in rels:
        if value.startswith(rel):
            return True
    return False


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
    stmt = Influence(subj_concept, obj_concept)
    stmt.obj_delta['polarity'] = pol
    return stmt
