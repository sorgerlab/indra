import unittest
import json

from db_rest_api import api
from datetime import datetime
from indra.statements import stmts_from_json
from itertools import combinations

TIMELIMIT = 1

class DbApiTestCase(unittest.TestCase):

    def setUp(self):
        api.app.testing = True
        self.app = api.app.test_client()

    def tearDown(self):
        pass

    def __time_get_query(self, query_str):
        start_time = datetime.now()
        resp = self.app.get('/statements/?%s' % query_str)
        t_delta = datetime.now() - start_time
        dt = t_delta.seconds + t_delta.microseconds/1e6
        print(dt)
        return resp, dt

    def __check_good_query(self, *args, **kwargs):
        query_str = '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()]
                             + list(args))
        resp, dt = self.__time_get_query(query_str)
        assert resp.status_code == 200, \
            'Got error code %d: \"%s\".' % (resp.status_code, resp.data.decode())
        json_stmts = json.loads(resp.data.decode('utf8'))
        assert len(json_stmts) is not 0, \
            'Did not get any statements.'
        stmts = stmts_from_json(json_stmts)
        assert len(stmts) == 1 or any([s.supports + s.supported_by for s in stmts]),\
            "Some statements lack support: %s." % str([str(s) for s in stmts if not s.supports + s.supported_by])
        assert all([not s1.matches(s2) for s1, s2 in combinations(stmts, 2)]),\
            "Some statements match: %s." % str([(s1, s2) for s1, s2 in combinations(stmts, 2) if s1.matches(s2)])
        assert dt <= TIMELIMIT, \
            "Query took %f seconds. Must be less than %f seconds." % (dt, TIMELIMIT)
        return resp

    def test_blank_response(self):
        """Test the response to an empty request."""
        resp, dt = self.__time_get_query('')
        assert resp.status_code == 400, \
            ('Got unexpected response with code %d: %s.'
             % (resp.status_code, resp.data.decode()))
        assert dt <= TIMELIMIT, \
            "Query took %f seconds. Must be less than %f seconds." % (dt, TIMELIMIT)

    def test_specific_query(self):
        """Test whether we can get a "fully" specified statement."""
        self.__check_good_query(object='MAP2K1', subject='MAPK1', type='Phosphorylation')

    def test_query_with_two_agents(self):
        """Test a query were the roles of the agents are not given."""
        self.__check_good_query('agent=MAP2K1', 'agent=MAPK1', type='Phosphorylation')

    def test_query_with_other(self):
        """Test that we can get an ActiveForm."""
        self.__check_good_query(agent='MAPK1', type='ActiveForm')

    def test_bad_camel(self):
        self.__check_good_query(agent='MAPK1', type='acTivefOrm')
        

if __name__ == '__main__':
    unittest.main()

