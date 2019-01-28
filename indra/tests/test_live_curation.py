import os
import json
import unittest
from indra.statements import *
from indra.tools.live_curation import app, corpora


def _make_corpus():
    test_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                  'test_corpus.json')
    ev1 = Evidence(source_api='eidos', text='A',
                   epistemics={'found_by': 'ported_syntax_1_verb-Causal'})
    ev2 = Evidence(source_api='eidos', text='B',
                   epistemics={'found_by': 'dueToSyntax2-Causal'})
    ev3 = Evidence(source_api='hume', text='C')
    ev4 = Evidence(source_api='cwms', text='D')
    ev5 = Evidence(source_api='sofia', text='E')
    stmt1 = Influence(Concept('x'), Concept('y'), evidence=[ev1, ev2])
    stmt2 = Influence(Concept('x'), Concept('y'), evidence=[ev1, ev3])
    stmt3 = Influence(Concept('x'), Concept('y'), evidence=[ev3, ev4, ev5])
    stmt4 = Influence(Concept('x'), Concept('y'), evidence=[ev5])
    stmts_to_json_file([stmt1, stmt2, stmt3, stmt4], test_file_path)

corpora['1'] = _make_corpus()

class LiveCurationTestCase(unittest.TestCase):
    def setUp(self):
        _make_corpus()
        app.testing = True
        self.app = app.test_client()

    # Tests ==================
    def test_alive(self):
        resp = self.app.post('update_beliefs',
                             data=json.dumps({'corpus_id': '1'}),
                             headers={'Content-Type': 'application/json'})
        assert resp.status_code == 200, resp
