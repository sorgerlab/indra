import json
import requests
from os import path
from datetime import datetime

from nose.plugins.attrib import attr

from indra.statements import stmts_from_json, Phosphorylation

BASE = 'http://api.indra.bio:8000/'
HERE = path.dirname(path.abspath(__file__))


def _call_api(method, route, *args, **kwargs):
    route = route.lstrip('/')
    req_meth = getattr(requests, method)
    start = datetime.now()
    print("Submitting request to '%s' at %s." % ('/' + route, start))
    print("\targs:", args)
    print("\tkwargs:", kwargs)
    res = req_meth(BASE + route, *args, **kwargs)
    end = datetime.now()
    print("Got result with %s at %s after %s seconds."
          % (res.status_code, end, (end-start).total_seconds()))
    assert res.status_code == 200, res.status_code
    return res


@attr('webservice')
def test_responsive():
    res = _call_api('get', '')
    assert res.content.startswith(b'This is the INDRA REST API.'), \
        "Unexpected content: %s" % res.content


@attr('webservice')
def test_options():
    res = _call_api('options', '')
    assert res.content == b'{}', \
        "Unexpected content: %s" % res.content


@attr('webservice', 'notravis')
def test_trips_process_text():
    res = _call_api('post', 'trips/process_text',
                    json={'text': 'MEK phosphorylates ERK.'})
    res_json = res.json()
    assert 'statements' in res_json.keys(), res_json
    stmts = stmts_from_json(res_json['statements'])
    assert len(stmts) == 1, len(stmts)
    stmt = stmts[0]
    assert isinstance(stmt, Phosphorylation), type(stmt)
    assert stmt.enz.name == 'MEK', stmt.enz
    assert stmt.sub.name == 'ERK', stmt.sub


@attr('webservice')
def test_trips_process_xml():
    test_ekb_path = path.join(HERE, 'trips_ekbs',
                              'MEK_increases_the_phosphorylation_of_ERK.ekb')
    with open(test_ekb_path, 'r') as f:
        xml_str = f.read()
    res = _call_api('post', 'trips/process_xml', json={'xml_str': xml_str})
    res_json = res.json()
    assert 'statements' in res_json.keys(), res_json
    stmts = stmts_from_json(res_json['statements'])
    assert len(stmts) == 1, len(stmts)
    stmt = stmts[0]
    assert isinstance(stmt, Phosphorylation), type(stmt)
    assert stmt.enz.name == 'MEK', stmt.enz
    assert stmt.sub.name == 'ERK', stmt.sub


@attr('webservice', 'notravis')
def test_reach_process_text():
    res = _call_api('post', 'reach/process_text',
                    json={'text': 'MEK phosphorylates ERK.'})
    res_json = res.json()
    assert 'statements' in res_json.keys(), res_json
    stmts = stmts_from_json(res_json['statements'])
    assert len(stmts) == 1, len(stmts)
    stmt = stmts[0]
    assert isinstance(stmt, Phosphorylation), type(stmt)
    assert stmt.enz.name == 'MEK', stmt.enz
    assert stmt.sub.name == 'ERK', stmt.sub


@attr('webservice')
def test_reach_process_json():
    # TODO: Add test of reach process json
    return


STMT_JSON = {
    "id": "acc6d47c-f622-41a4-8ae9-d7b0f3d24a2f",
    "type": "Complex",
    "members": [
        {"db_refs": {"TEXT": "MEK", "FPLX": "MEK"}, "name": "MEK"},
        {"db_refs": {"TEXT": "ERK", "FPLX": "ERK"}, "name": "ERK"}
    ],
    "sbo": "http://identifiers.org/sbo/SBO:0000526",
    "evidence": [{"text": "MEK binds ERK", "source_api": "trips"}]
}


@attr('webservice')
def test_assemblers_cyjs():
    stmt_str = json.dumps({'statements': [STMT_JSON]})
    res = _call_api('post', 'assemblers/cyjs', stmt_str)
    res_json = res.json()
    assert len(res_json['edges']) == 1, len(res_json['edges'])
    assert len(res_json['nodes']) == 2, len(res_json['nodes'])
    return


@attr('webservice')
def test_assemblers_pysb_no_format():
    stmt_str = json.dumps({'statements': [STMT_JSON]})
    res = _call_api('post', 'assemblers/pysb', stmt_str)
    res_json = res.json()
    assert 'model' in res_json.keys()
    assert res_json['model'] is not None
    assert 'MEK' in res_json['model'], res_json['model']
    return


@attr('webservice')
def test_assemblers_pysb_kappa_img_format():
    for exp_format in ['kappa_im', 'kappa_cm']:
        print("Testing", exp_format)
        stmt_str = json.dumps({'statements': [STMT_JSON],
                               'export_format': exp_format})
        res = _call_api('post', 'assemblers/pysb', stmt_str)
        res_json = res.json()
        assert 'image' in res_json.keys()
        assert 'image' in res_json.keys()
        assert res_json['image'] is not None
    return


@attr('webservice')
def test_assemblers_pysb_kappa_other_formats():
    # All the formats defined in PysbAssembler.export_model doc string.
    formats = ['bngl', 'kappa', 'sbml', 'matlab', 'mathematica',
               'potterswheel']
    for exp_format in formats:
        print("Testing", exp_format)
        stmt_str = json.dumps({'statements': [STMT_JSON],
                               'export_format': exp_format})
        res = _call_api('post', 'assemblers/pysb', stmt_str)
        res_json = res.json()
        assert 'model' in res_json.keys()
        assert res_json['model'] is not None
        assert 'MEK' in res_json['model'], res_json['model']
    return


@attr('webservice')
def test_assemblers_cx():
    stmt_str = json.dumps({'statements': [STMT_JSON]})
    res = _call_api('post', 'assemblers/cx', stmt_str)
    res_json = res.json()
    assert 'model' in res_json.keys()
    assert res_json['model'] is not None
    assert 'MEK' in res_json['model'], res_json['model']


@attr('webservice')
def test_assemblers_graph():
    stmt_str = json.dumps({'statements': [STMT_JSON]})
    res = _call_api('post', 'assemblers/graph', stmt_str)
    res_json = res.json()
    assert 'model' in res_json.keys()
    assert res_json['model'] is not None
    assert 'MEK' in res_json['model'], res_json['model']


@attr('webservice')
def test_assemblers_english():
    stmt_str = json.dumps({'statements': [STMT_JSON]})
    res = _call_api('post', 'assemblers/english', stmt_str)
    res_json = res.json()
    assert 'sentences' in res_json.keys()
    assert len(res_json['sentences']) == 1, len(res_json['sentences'])
    sentence = list(res_json['sentences'].values())[0]
    assert 'MEK' in sentence, sentence


@attr('webservice')
def test_assemblers_loopy():
    stmt_jsons = [{
            "id": "acc6d47c-f622-41a4-8ae9-d7b0f3d24a2f",
            "type": "Phosphorylation",
            "enz": {"db_refs": {"TEXT": "MEK", "FPLX": "MEK"}, "name": "MEK"},
            "sub": {"db_refs": {"TEXT": "ERK", "FPLX": "ERK"}, "name": "ERK"},
            "sbo": "http://identifiers.org/sbo/SBO:0000526",
            "evidence": [{"text": "MEK phosphorylates ERK", "source_api": "trips"}]
        },
        {
            "id": "bcc6d47c-f622-41a4-8ae9-d7b0f3d24a2f",
            "type": "Activation",
            "subj": {"db_refs": {"TEXT": "ERK", "FPLX": "ERK"}, "name": "ERK"},
            "obj": {"db_refs": {"TEXT": "EGFR", "HGNC": "3236"}, "name": "EGFR"},
            "sbo": "http://identifiers.org/sbo/SBO:0000526",
            "evidence": [{"text": "ERK activates EGFR", "source_api": "trips"}]
        }
    ]
    stmt_str = json.dumps({'statements': stmt_jsons})
    res = _call_api('post', 'assemblers/sif/loopy', stmt_str)
    res_json = res.json()
    assert 'loopy_url' in res_json.keys()
    assert "ERK" in res_json['loopy_url']


@attr('webservice')
def test_pipeline():
    p = [{'function': 'filter_grounded_only'},
         {'function': 'run_preassembly',
          'kwargs': {'return_toplevel': False}}]
    json_str = json.dumps({'statements': [STMT_JSON], 'pipeline': p})
    res = _call_api('post', 'preassembly/pipeline', json_str)
    res_json = res.json()
    assert 'statements' in res_json
    assert len(res_json['statements']) == 1
