import json
import unittest
from indra.statements import *
from indra.tools.live_curation import app, corpora, Corpus, scorer, \
    default_priors, wm_scorer


def _make_corpus():
    ev1 = Evidence(source_api='eidos', text='A',
                   annotations={'found_by': 'ported_syntax_1_verb-Causal'})
    ev2 = Evidence(source_api='eidos', text='B',
                   annotations={'found_by': 'dueToSyntax2-Causal'})
    ev3 = Evidence(source_api='hume', text='C')
    ev4 = Evidence(source_api='cwms', text='D')
    ev5 = Evidence(source_api='sofia', text='E')
    ev6 = Evidence(source_api='sofia', text='F')
    stmt1 = Influence(Concept('x'), Concept('y'), evidence=[ev1, ev2])
    stmt2 = Influence(Concept('x'), Concept('y'), evidence=[ev1, ev3])
    stmt3 = Influence(Concept('x'), Concept('y'), evidence=[ev3, ev4, ev5])
    stmt4 = Influence(Concept('x'), Concept('y'), evidence=[ev5])
    stmt5 = Influence(Concept('x'), Concept('y'), evidence=[ev6])
    stmt1.uuid = '1'
    stmt2.uuid = '2'
    stmt3.uuid = '3'
    stmt4.uuid = '4'
    stmt5.uuid = '5'
    return Corpus([stmt1, stmt2, stmt3, stmt4, stmt5])


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
        self._reset()
        resp = self._send_request({'corpus_id': '1'})
        assert resp.status_code == 200, resp

    def test_bad_corpus(self):
        self._reset()
        resp = self._send_request({'corpus_id': '2'})
        assert resp.status_code == 400, resp

    def test_no_curation(self):
        self._reset()
        resp = self._send_request({'corpus_id': '1',
                                   'return_beliefs': True})
        res = json.loads(resp.data.decode('utf-8'))
        expected = {'1': 0.9167547741034001,
                    '2': 0.8968421052631579,
                    '3': 0.957125,
                    '4': 0.65,
                    '5': 0.65}
        assert close_enough(res, expected), (res, expected)

    def test_eid_rule1_incorrect(self):
        self._reset()
        resp = self._send_request({'corpus_id': '1',
                                   'curations': {'1': 0},
                                   'return_beliefs': True})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        expected = {'1': 0,
                    '2': 0.8917525773195876,
                    '3': 0.957125,
                    '4': 0.65,
                    '5': 0.65}
        assert close_enough(res, expected), (res, expected)

    def test_eid_rule1_incorrect_again(self):
        self._reset()
        resp = self._send_request({'corpus_id': '1',
                                   'curations': {'1': 0},
                                   'return_beliefs': True})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        expected = {'1': 0,
                    '2': 0.8917525773195876,
                    '3': 0.957125,
                    '4': 0.65,
                    '5': 0.65}
        assert close_enough(res, expected), (res, expected)


    def test_eid_rule1_correct(self):
        self._reset()
        resp = self._send_request({'corpus_id': '1',
                                   'curations': {'1': 1},
                                   'return_beliefs': True})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        expected = {'1': 1,
                    '2': 0.8979166666666667,
                    '3': 0.957125,
                    '4': 0.65,
                    '5': 0.65}
        assert close_enough(res, expected), (res, expected)


    def test_eid_rule2_correct(self):
        self._reset()
        resp = self._send_request({'corpus_id': '1',
                                   'curations': {'2': 1},
                                   'return_beliefs': True})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert res == {'1': 0.9171718289085546,
                       '2': 1,
                       '3': 0.9591666666666667,
                       '4': 0.65,
                       '5': 0.65}, res

    def test_hume_correct(self):
        self._reset()
        resp = self._send_request({'corpus_id': '1',
                                   'curations': {'3': 0},
                                   'return_beliefs': True})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert close_enough(res, {'1': 0.9167547741034001,
                                  '2': 0.887719298245614,
                                  '3': 0,
                                  '4': 0.6190476190476191,
                                  '5': 0.6190476190476191}), res

    def test_sofia_incorrect(self):
        self._reset()
        resp = self._send_request({'corpus_id': '1',
                                   'curations': {'4': 0},
                                   'return_beliefs': True})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert res == {'1': 0.9167547741034001,
                       '2': 0.8968421052631579,
                       '3': 0.9533333333333334,
                       '4': 0,
                       '5': 0.6190476190476191}, res
        resp = self._send_request({'corpus_id': '1',
                                   'curations': {'5': 0},
                                   'return_beliefs': True})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert close_enough(res, {'1': 0.9167547741034001,
                                  '2': 0.8968421052631579,
                                  '3': 0.9498863636363637,
                                  '4': 0,
                                  '5': 0}), res

    def _reset(self):
        global corpora, scorer
        corpora['1'] = _make_corpus()
        scorer = wm_scorer.get_eidos_bayesian_scorer(default_priors)
        self.app = app.test_client()


def close_enough(probs, ref):
    for k, v in probs.items():
        if abs(ref[k] - probs[k]) > 0.001:
            return False
    return True
