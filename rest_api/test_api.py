import os
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
    assert len(stmts) == 1


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


def test_filter_grounded_only_score():
    db_refs1 = {'UN': [('x', 0.1)]}
    db_refs2 = {'UN': [('x', 0.5)]}
    st1 = Influence(Concept('a', db_refs=db_refs1),
                    Concept('b', db_refs=db_refs1))
    st2 = Influence(Concept('a', db_refs=db_refs1),
                    Concept('b', db_refs=db_refs2))
    st3 = Influence(Concept('a', db_refs=db_refs2),
                    Concept('b', db_refs=db_refs2))
    stmts_json = stmts_to_json([st1, st2, st3])
    url = base_url + '/preassembly/filter_grounded_only'
    res = requests.post(url, json={'statements': stmts_json,
                                   'score_threshold': 0.3})
    res_json = res.json()
    stmts_json = res_json.get('statements')
    assert len(stmts_json) == 1, len(stmts_json)


def test_loopy():
    url = base_url + '/reach/process_text'
    res = requests.post(url, json={'text': 'MEK activates ERK.'})
    url = base_url + '/assemblers/sif/loopy'
    res = requests.post(url, json=res.json())
    res_json = res.json()
    print(res_json.get('loopy_url'))


def test_cwms():
    url = base_url + '/cwms/process_text'
    res = requests.post(url, json={'text': 'Hunger causes displacement.'})
    res_json = res.json()
    stmts_json = res_json.get('statements')
    stmts = stmts_from_json(stmts_json)
    assert len(stmts) == 1


def test_hume():
    from indra.tests.test_hume import test_file_new_simple
    with open(test_file_new_simple, 'r') as fh:
        test_jsonld = fh.read()
    url = base_url + '/hume/process_jsonld'
    res = requests.post(url, json={'jsonld': test_jsonld})
    print(res.content)
    res_json = res.json()
    stmts_json = res_json.get('statements')
    stmts = stmts_from_json(stmts_json)
    assert len(stmts) == 1


def test_eidos_json():
    from indra.tests.test_eidos import test_jsonld
    with open(test_jsonld, 'r') as fh:
        test_json = fh.read()
    url = base_url + '/eidos/process_json'
    res = requests.post(url, json={'json': test_json})
    print(res.content)
    res_json = res.json()
    stmts_json = res_json.get('statements')
    stmts = stmts_from_json(stmts_json)
    assert len(stmts) == 1


def test_belief_filter():
    st1 = Influence(Concept('a'), Concept('b'))
    st1.belief = 0.2
    st2 = Influence(Concept('a'), Concept('b'))
    st2.belief = 0.5
    st3 = Influence(Concept('a'), Concept('b'))
    st3.belief = 0.8
    stmts_json = stmts_to_json([st1, st2, st3])
    url = base_url + '/preassembly/filter_belief'
    res = requests.post(url, json={'statements': stmts_json,
                                   'belief_cutoff': 0.6})
    res_json = res.json()
    stmts_json = res_json.get('statements')
    assert len(stmts_json) == 1, len(stmts_json)
