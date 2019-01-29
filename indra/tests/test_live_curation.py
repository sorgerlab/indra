import json
import unittest
from indra.statements import *
from indra.tools.live_curation import app, corpora, Corpus


def _make_corpus():
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
    stmt1.uuid = '1'
    stmt2.uuid = '2'
    stmt3.uuid = '3'
    stmt4.uuid = '4'
    return Corpus([stmt1, stmt2, stmt3, stmt4])


corpora['1'] = _make_corpus()


class LiveCurationTestCase(unittest.TestCase):
    def setUp(self):
        _make_corpus()
        app.testing = True
        self.app = app.test_client()

    def _send_request(self, req_dict):
        resp = self.app.post('update_beliefs',
                             data=json.dumps(req_dict),
                             headers={'Content-Type': 'application/json'})
        return resp

    # Tests ==================
    def test_alive(self):
        resp = self._send_request({'corpus_id': '1'})
        assert resp.status_code == 200, resp

    def test_bad_corpus(self):
        resp = self._send_request({'corpus_id': '2'})
        assert resp.status_code == 400, resp

    def test_eidos_rule1(self):
        resp = self._send_request({'corpus_id': '1',
                                   'curations': {'1': 0},
                                   'return_beliefs': True})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert res['1'] == 0
