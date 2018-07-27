import os
from indra.sources import eidos
from indra.statements import Influence
import requests
import json
from indra.assemblers.cag_assembler import CAGAssembler
from indra.assemblers.cx_assembler import CxAssembler
from indra.assemblers.pysb_assembler import PysbAssembler


path_this = os.path.dirname(os.path.abspath(__file__))
test_jsonld = os.path.join(path_this, 'eidos_test.jsonld')


def __get_remote_jsonld():
    res = requests.get('https://raw.githubusercontent.com/clulab/eidos/master/'
                       'example_output/example-0.2.2.jsonld')
    assert res.status_code is 200, "Could not get example json from remote."
    example_json = json.loads(res.content.decode('utf-8'))
    return example_json


def __get_stmts_from_remote_jsonld():
    ex_json = __get_remote_jsonld()
    ep = eidos.process_json_ld(ex_json)
    assert ep is not None, 'Failed to handle json with eidos processor.'
    assert len(ep.statements), 'Did not get statements from json.'
    return ep.statements


def test_process_text():
    ep = eidos.process_text('The cost of fuel decreases water trucking.',
                            out_format='json_ld')
    assert ep is not None
    assert len(ep.statements) == 1
    stmt = ep.statements[0]
    assert isinstance(stmt, Influence)
    assert stmt.subj.name == 'cost fuel'
    assert stmt.obj.name == 'water trucking'
    assert stmt.obj_delta.get('polarity') == -1
    assert(stmt.evidence[0].annotations['found_by']
           == 'ported_syntax_1_verb-Causal')


def test_process_text_json_ld():
    ep = eidos.process_text('The cost of fuel decreases water trucking.',
                            out_format='json_ld')
    assert ep is not None
    assert len(ep.statements) == 1
    stmt = ep.statements[0]
    assert isinstance(stmt, Influence)
    assert stmt.subj.name == 'cost fuel'
    assert stmt.obj.name == 'water trucking'
    assert stmt.obj_delta.get('polarity') == -1
    assert(stmt.evidence[0].annotations['found_by']
           == 'ported_syntax_1_verb-Causal')
    assert 'TEXT' in stmt.subj.db_refs
    assert 'TEXT' in stmt.obj.db_refs
    # assert 'UN' in stmt.subj.db_refs
    # assert 'UN' in stmt.obj.db_refs
    # FIXME: once groundings are propagated well from offline reading
    # this should work
    # assert len(stmt.subj.db_refs['UN']) > 5
    # assert len(stmt.obj.db_refs['UN']) > 5
    # Make sure sanitization works
    sanitized = ep._sanitize('-LRB-something-RRB-')
    assert sanitized == '(something)'


def test_process_json_ld_file():
    ep = eidos.process_json_ld_file(test_jsonld)
    assert len(ep.statements) == 1
    assert 'UN' in ep.statements[0].subj.db_refs
    assert 'UN' in ep.statements[0].obj.db_refs

def test_eidos_to_cag():
    stmts = __get_stmts_from_remote_jsonld()
    ca = CAGAssembler()

    # Make sure these don't error
    ca.add_statements(stmts)
    ca.make_model()
    ca.export_to_cytoscapejs()
    return


def test_eidos_to_cx():
    stmts = __get_stmts_from_remote_jsonld()
    cx = CxAssembler()

    # Make sure these don't error
    cx.add_statements(stmts)
    cx.make_model()
    test_fname = 'test_cag_to_cx.cx'
    try:
        cx.save_model(test_fname)
        assert os.path.exists(test_fname), "Failed to create cx file."
    finally:
        if os.path.exists(test_fname):
            os.remove(test_fname)
    return


def test_eidos_to_pysb():
    stmts = __get_stmts_from_remote_jsonld()
    pa = PysbAssembler()

    # Make sure these don't error
    pa.add_statements(stmts)
    pa.make_model()
    for fmt in ['kappa', 'sbml', 'sbgn']:
        exp_str = pa.export_model(fmt)
        assert exp_str, "Got no exported model from eidos->psyb to %s." % fmt
    return
