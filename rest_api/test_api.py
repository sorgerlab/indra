import requests
from indra.statements import *

base_url = 'http://localhost:8080'

def test_filter_by_type():
    st1 = Phosphorylation(Agent('a'), Agent('b'))
    st2 = Complex([Agent('a'), Agent('b')])
    stmts_json = stmts_to_json([st1, st2])
    url = base_url + '/preassembly/filter_by_type'
    data = {'statements': stmts_json,
            'type': 'phosphorylation'}
    res = requests.post(url, json=data)
    res_json = res.json()
    stmts_json = res_json.get('statements')
    stmts = stmts_from_json(stmts_json)
    assert(len(stmts) == 1)

def test_filter_grounded_only():
    a = Agent('a', db_refs={'HGNC': '1234'})
    b = Agent('b', db_refs={'HGNC': '1235'})
    c = Agent('c', db_refs={'TEXT': 'c'})
    d = Agent('d', db_refs={})
    st1 = Phosphorylation(a, b)
    st2 = Phosphorylation(a, c)
    st3 = Phosphorylation(a, d)
    stmts_json = stmts_to_json([st1, st2, st3])
    url = base_url + '/preassembly/filter_grounded_only'
    data = {'statements': stmts_json,
            'type': 'phosphorylation'}
    res = requests.post(url, json=data)
    res_json = res.json()
    stmts_json = res_json.get('statements')
    stmts = stmts_from_json(stmts_json)
    assert(len(stmts) == 1)

def test_loopy():
    url = base_url + '/reach/process_text'
    res = requests.post(url, json={'text': 'MEK activates ERK.'})
    url = base_url + '/assemblers/sif/loopy'
    res = requests.post(url, json=res.json())
    res_json = res.json()
    print(res_json.get('loopy_url'))
