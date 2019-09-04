import openpyxl
from indra.statements import *


def process_sheet(sheet):
    stmts = []
    header = [cell.value for cell in next(sheet.rows)]
    for row in sheet.rows:
        row_values = [r.value for r in row]
        if not any(row_values):
            break
        row_dict = {h: v for h, v in zip(header, row_values)}
        stmt = make_stmt(row_dict)
        if stmt:
            stmts.append(stmt)
    return stmts


def make_stmt(row_dict):
    # Make subject
    subj_concept = Concept(
        name=row_dict['Source/Cause (Factor A)'],
        db_refs=_make_grounding(
            row_dict['Source/Cause (Factor A)'],
            row_dict['Source/Cause node (WM ontology node)']))
    # Handle case where cause is missing
    if subj_concept.name is None:
        return None
    if row_dict['Source/Cause polarity']:
        if row_dict['Source/Cause polarity'].lower() == 'increase':
            subj_pol = 1
        elif row_dict['Source/Cause polarity'].lower() == 'decrease':
            subj_pol = -1
    else:
        subj_pol = None
    subj_time = TimeContext(
        text=str(row_dict['Original temporal text for cause']))
    subj_context = WorldContext(time=subj_time, geo_location=RefContext(
        name=row_dict['CauseLocation']))
    subj = Event(subj_concept, delta=QualitativeDelta(polarity=subj_pol),
                 context=subj_context)

    # Make object
    obj_concept = Concept(
        row_dict['Target/Effect (Factor B)'],
        db_refs=_make_grounding(row_dict['Target/Effect (Factor B)'],
                                row_dict['Target/Effect (WM ontology node)']))
    # Handle case where effect is missing
    if subj_concept.name is None:
        return None
    obj_time = TimeContext(
        text=str(row_dict['Original temporal text for effect']))
    if row_dict['Target/Effect polarity']:
        if row_dict['Target/Effect polarity'].lower() == 'increase':
            obj_pol = 1
        elif row_dict['Target/Effect polarity'].lower() == 'decrease':
            obj_pol = -1
    else:
        obj_pol = None
    # If both source and destination locations are present, make Migration obj
    if row_dict['Effect SourceLocation'] and \
            row_dict['Effect DestinationLocation']:
        obj_context = MovementContext(locations=[
            {'location': RefContext(name=row_dict['Effect SourceLocation']),
             'role': 'origin'},
            {'location': RefContext(name=row_dict['Effect DestinationLocation']),
             'role': 'destination'}], time=obj_time)
        obj = Migration(
            obj_concept, delta=QuantitativeState(
                entity='person', value=row_dict['Relation strength estimate'],
                unit=row_dict['Relation strength unit'], polarity=obj_pol),
            context=obj_context)
    # If only one or no location is present, make Event object
    else:
        if row_dict['Effect SourceLocation']:
            obj_loc = RefContext(name=row_dict['Effect SourceLocation'])
        elif row_dict['Effect DestinationLocation']:
            obj_loc = RefContext(name=row_dict['Effect DestinationLocation'])
        else:
            obj_loc = None
        obj_context = WorldContext(time=obj_time, geo_location=obj_loc)
        obj = Event(obj_concept, delta=QualitativeDelta(polarity=obj_pol),
                    context=obj_context)
    # Make evidence
    anns = {'provenance': row_dict['Provenance (file id)'],
            'data_format': row_dict['Data format'],
            'group_id': row_dict['Group ID']}
    if anns['data_format'] == 'text':
        text = row_dict['Sentence(s)/Figure/Table']
    else:
        anns['table_figure_name'] = row_dict['Sentence(s)/Figure/Table']
        text = None
    ev = Evidence(text=text, annotations=anns, source_api='assertion')
    stmt = Influence(subj, obj, [ev])
    return stmt


def _make_grounding(text, ont_str):
    # Make a list of tuples of scored ontology concepts from a single string
    grounding = {'TEXT': text}
    if ont_str:
        wm_grounding = []
        ont_list = ont_str[1:-1].split(') (')
        for entry in ont_list:
            ont_concept, score = entry.split(',')
            wm_grounding.append((ont_concept, float(score)))
        grounding['WM'] = wm_grounding
    return grounding


def _make_event_migration(concept):
    grounding = concept.get_grounding()
    if grounding:
        if grounding[1].startswith(
                'wm/concept/causal_factor/social_and_political/migration'):
            return Migration(concept)
        else:
            return Event(concept)
    return Event(concept)


def process_workbook(fname):
    wb = openpyxl.load_workbook(fname, read_only=True)
    sheets = wb.sheetnames
    cag_sheets = [s for s in sheets if 'CAG' in s]

    statements = []
    for sheet_name in cag_sheets:
        sheet = wb[sheet_name]
        new_stmts = process_sheet(sheet)
        statements += new_stmts
    return statements
