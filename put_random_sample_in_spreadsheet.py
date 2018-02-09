"""Takes the pickled random sample and outputs a spreadsheet for
curation"""

from indra.statements import *

import pickle

fname = 'reach_random_sample_size300.p'
sample = pickle.load( open( fname, "rb" ) )

for statement in sample:
    statement_str = repr(statement)

    if isinstance(statement, Translocation):
        subj = statement.agent
        subj_str = str(subj)
        subj_text = subj.db_refs['TEXT']

        obj_str = str(statement.from_location) + '->' + str(statement.to_location)
        obj_text = ''
    elif isinstance(statement, Modification):
        subj = statement.enz
        obj = statement.sub

        subj_str = str(subj)
        try:
            subj_text = subj.db_refs['TEXT']
        except:
            subj_text = ''
        obj_str = str(obj)
        try:
            obj_text = obj.db_refs['TEXT']
        except:
            obj_text = ''
    elif isinstance(statement, Complex):
        members = statement.members

        if len(members) > 2:
            obj_str = ''
            obj_text = ''
            subj_str = ''
            subj_text = ''
        else:
            subj = members[0]

            subj_str = str(subj)
            subj_text = subj.db_refs['TEXT']

            if len(members) == 2:
                obj = members[1]
                obj_str = str(obj)
                obj_text = obj.db_refs['TEXT']
            else:
                obj_str = ''
                obj_text = ''

    else:
        subj = statement.subj
        subj_str = str(subj)
        subj_text = subj.db_refs['TEXT']

        obj = statement.obj
        obj_str = str(obj)
        obj_text = obj.db_refs['TEXT']

    text = statement.evidence[0].text
    reach_rule = statement.evidence[0].annotations['found_by']

    row = [statement_str, obj_str, obj_text, subj_str, subj_text, text, reach_rule]
    print('~'.join(row))
    print(statement.evidence[0].pmid)

