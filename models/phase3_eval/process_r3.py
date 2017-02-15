from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
from indra.statements import Statement, Evidence

active_forms_file = 'sources/cure-active-forms-2017-01-30.txt'

def read_jsons(fname):
    all_json = []
    with open(fname, 'rt') as fh:
        lin = fh.readlines()
    lin = [l.strip() for l in lin if l.strip()]
    for l in lin:
        lin_json = json.loads(l)
        all_json.append(lin_json[0])
    return all_json

def read_stmts(fname):
    jsons = read_jsons(fname)
    stmts = []
    for js in jsons:
        st = Statement.from_json(json.dumps(js))
        stmts.append(st)
    for st in stmts:
        # Correct evidence
        if isinstance(st.evidence, Evidence):
            st.evidence = [st.evidence]
        for ev in st.evidence:
            ev.source_api = 'r3'
            ev.pmid = None
            ev.source_id = None
            ev.annotations = {}
            ev.epistemics = {}
        # Correct supports/supported by
        st.supports = []
        st.supported_by = []
        st.belief = 1
    return stmts

if __name__ == '__main__':
    stmts = read_stmts(active_forms_file)
