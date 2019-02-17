from indra.statements import *
from indra.tools.live_curation import Corpus, LiveCurator


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


def test_no_curation():
    curator = LiveCurator(corpora={'1': _make_corpus()})
    curator.submit_curation(corpus_id='1', curations={})
    beliefs = curator.update_beliefs(corpus_id='1')
    expected = {'1': 0.91675,
                '2': 0.8968,
                '3': 0.957125,
                '4': 0.65,
                '5': 0.65}
    assert close_enough(beliefs, expected), (beliefs, expected)


def test_eid_rule1_incorrect():
    curator = LiveCurator(corpora={'1': _make_corpus()})
    curator.submit_curation(corpus_id='1', curations={'1': 0})
    expected = {'1': 0,
                '2': 0.8942,
                '3': 0.957125,
                '4': 0.65,
                '5': 0.65}
    beliefs = curator.update_beliefs(corpus_id='1')
    assert close_enough(beliefs, expected), (beliefs, expected)

    # Submit another curation
    curator.submit_curation(corpus_id='1', curations={'1': 0})
    expected = {'1': 0,
                '2': 0.8917,
                '3': 0.957125,
                '4': 0.65,
                '5': 0.65}
    beliefs = curator.update_beliefs(corpus_id='1')
    assert close_enough(beliefs, expected), (beliefs, expected)


def test_eid_rule1_correct():
    curator = LiveCurator(corpora={'1': _make_corpus()})
    curator.submit_curation(corpus_id='1', curations={'1': 1})
    expected = {'1': 1,
                '2': 0.8979,
                '3': 0.957125,
                '4': 0.65,
                '5': 0.65}
    beliefs = curator.update_beliefs(corpus_id='1')
    assert close_enough(beliefs, expected), (beliefs, expected)


def test_eid_rule2_correct():
    curator = LiveCurator(corpora={'1': _make_corpus()})
    curator.submit_curation(corpus_id='1', curations={'2': 1})
    expected = {'1': 0.91717,
                '2': 1,
                '3': 0.95916,
                '4': 0.65,
                '5': 0.65}
    beliefs = curator.update_beliefs(corpus_id='1')
    assert close_enough(beliefs, expected), (beliefs, expected)


def test_hume_incorrect():
    curator = LiveCurator(corpora={'1': _make_corpus()})
    curator.submit_curation(corpus_id='1', curations={'3': 0})
    expected = {'1': 0.91675,
                '2': 0.88772,
                '3': 0,
                '4': 0.61904,
                '5': 0.61904}
    beliefs = curator.update_beliefs(corpus_id='1')
    assert close_enough(beliefs, expected), (beliefs, expected)


def test_sofia_incorrect():
    curator = LiveCurator(corpora={'1': _make_corpus()})
    curator.submit_curation(corpus_id='1', curations={'4': 0})
    expected = {'1': 0.91675,
                '2': 0.89684,
                '3': 0.95333,
                '4': 0.0,
                '5': 0.61904}
    beliefs = curator.update_beliefs(corpus_id='1')
    assert close_enough(beliefs, expected), (beliefs, expected)

    curator.submit_curation(corpus_id='1', curations={'5': 0})
    expected = {'1': 0.91675,
                '2': 0.89684,
                '3': 0.94988,
                '4': 0,
                '5': 0}
    beliefs = curator.update_beliefs(corpus_id='1')
    assert close_enough(beliefs, expected), (beliefs, expected)


"""
FIXME: this Flask-app test suite is not working
class LiveCurationTestCase(unittest.TestCase):
    def setUp(self):
        _make_corpus()
        app.testing = True
        self.app = app.test_client()

    def _send_request(self, endpoint, req_dict):
        resp = self.app.post(endpoint,
                             data=json.dumps(req_dict),
                             headers={'Content-Type': 'application/json'})
        return resp

    # Tests ==================
    def test_alive(self):
        resp = self._send_request('submit_curation', {'corpus_id': '1'})
        assert resp.status_code == 200, resp

    def test_bad_corpus(self):
        resp = self._send_request('submit_curation', {'corpus_id': '2'})
        assert resp.status_code == 400, resp

    def test_no_curation(self):
        self._send_request('submit_curation', {'corpus_id': '1'})
        resp = self._send_request('update_beliefs', {'corpus_id': '1'})
        res = json.loads(resp.data.decode('utf-8'))
        expected = {'1': 0.9167547741034001,
                    '2': 0.8968421052631579,
                    '3': 0.957125,
                    '4': 0.65,
                    '5': 0.65}
        assert close_enough(res, expected), (res, expected)

    def test_eid_rule1_incorrect(self):
        self._send_request('submit_curation', {'corpus_id': '1',
                                               'curations': {'1': 0}})
        resp = self._send_request('update_beliefs', {'corpus_id': '1'})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        expected = {'1': 0,
                    '2': 0.8917525773195876,
                    '3': 0.957125,
                    '4': 0.65,
                    '5': 0.65}
        assert close_enough(res, expected), (res, expected)

    def test_eid_rule1_incorrect_again(self):
        self._send_request('submit_curation', {'corpus_id': '1',
                                               'curations': {'1': 0}})
        resp = self._send_request('update_beliefs', {'corpus_id': '1'})
        assert resp.status_code == 200, resp
        res = json.loads(resp.data.decode('utf-8'))
        expected = {'1': 0,
                    '2': 0.8917525773195876,
                    '3': 0.957125,
                    '4': 0.65,
                    '5': 0.65}
        assert close_enough(res, expected), (res, expected)

    def test_eid_rule1_correct(self):
        resp = self._send_request('submit_curation', {'corpus_id': '1',
                                   'curations': {'1': 1}})
        resp = self._send_request('update_beliefs', {'corpus_id': '1'})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        expected = {'1': 1,
                    '2': 0.8979166666666667,
                    '3': 0.957125,
                    '4': 0.65,
                    '5': 0.65}
        assert close_enough(res, expected), (res, expected)

    def test_eid_rule2_correct(self):
        resp = self._send_request('submit_curation', {'corpus_id': '1',
                                   'curations': {'2': 1}})
        resp = self._send_request('update_beliefs', {'corpus_id': '1'})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert res == {'1': 0.9171718289085546,
                       '2': 1,
                       '3': 0.9591666666666667,
                       '4': 0.65,
                       '5': 0.65}, res

    def test_hume_correct(self):
        resp = self._send_request('submit_curation', {'corpus_id': '1',
                                   'curations': {'3': 0}})
        resp = self._send_request('update_beliefs', {'corpus_id': '1'})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert close_enough(res, {'1': 0.9167547741034001,
                                  '2': 0.887719298245614,
                                  '3': 0,
                                  '4': 0.6190476190476191,
                                  '5': 0.6190476190476191}), res

    def test_sofia_incorrect(self):
        resp = self._send_request('submit_curation', {'corpus_id': '1',
                                   'curations': {'4': 0}})
        resp = self._send_request('update_beliefs', {'corpus_id': '1'})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert res == {'1': 0.9167547741034001,
                       '2': 0.8968421052631579,
                       '3': 0.9533333333333334,
                       '4': 0,
                       '5': 0.6190476190476191}, res
        resp = self._send_request('submit_curation', {'corpus_id': '1',
                                   'curations': {'5': 0}})
        resp = self._send_request('update_beliefs', {'corpus_id': '1'})
        assert resp.status_code == 200
        res = json.loads(resp.data.decode('utf-8'))
        assert close_enough(res, {'1': 0.9167547741034001,
                                  '2': 0.8968421052631579,
                                  '3': 0.9498863636363637,
                                  '4': 0,
                                  '5': 0}), res
"""


def close_enough(probs, ref):
    for k, v in probs.items():
        if abs(ref[k] - probs[k]) > 0.0001:
            return False
    return True
