from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
from indra.statements import Statement

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
    return stmts

if __name__ == '__main__':
    stmts = read_stmts(active_forms_file)
