import unittest
import json
import sys

from itertools import combinations
from datetime import datetime

from indra.statements import stmts_from_json

from db_rest_api import api


TIMELIMIT = 1
SIZELIMIT = 4e7


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
        size = int(resp.headers['Content-Length'])
        raw_size = sys.getsizeof(resp.data)
        print("Raw size: %f, Compressed size: %f." % (raw_size/1e6, size/1e6))
        return resp, dt, size

    def __check_good_query(self, *args, **kwargs):
        check_stmts = kwargs.pop('check_stmts', True)
        time_limit = kwargs.pop('time_limit', TIMELIMIT)
        query_str = '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()]
                             + list(args))
        resp, dt, size = self.__time_get_query(query_str)
        assert resp.status_code == 200, \
            ('Got error code %d: \"%s\".'
             % (resp.status_code, resp.data.decode()))
        json_stmts = json.loads(resp.data.decode('utf8'))
        assert len(json_stmts) is not 0, \
            'Did not get any statements.'
        assert size <= SIZELIMIT, \
            ("Query took up %f MB. Must be less than %f MB."
             % (size/1e6, SIZELIMIT/1e6))
        stmts = stmts_from_json(json_stmts)
        assert any([s.supports + s.supported_by for s in stmts]),\
            ("Some statements lack support: %s."
             % str([str(s) for s in stmts if not s.supports + s.supported_by]))
        if check_stmts:
            assert all([not s1.matches(s2)
                        for s1, s2 in combinations(stmts, 2)]),\
                ("Some statements match: %s."
                 % str([(s1, s2) for s1, s2 in combinations(stmts, 2)
                        if s1.matches(s2)]))
        assert dt <= time_limit, \
            ("Query took %f seconds. Must be less than %f seconds."
             % (dt, time_limit))
        return resp

    def test_blank_response(self):
        """Test the response to an empty request."""
        resp, dt, size = self.__time_get_query('')
        assert resp.status_code == 400, \
            ('Got unexpected response with code %d: %s.'
             % (resp.status_code, resp.data.decode()))
        assert dt <= TIMELIMIT, \
            ("Query took %f seconds. Must be less than %f seconds."
             % (dt, TIMELIMIT))
        assert size <= SIZELIMIT, \
            "Query took up %f MB. Must be less than %f MB." % (size/1e6,
                                                               SIZELIMIT/1e6)

    def test_specific_query(self):
        """Test whether we can get a "fully" specified statement."""
        self.__check_good_query(object='MAP2K1', subject='MAPK1',
                                type='Phosphorylation')

    def test_query_with_two_agents(self):
        """Test a query were the roles of the agents are not given."""
        self.__check_good_query('agent=MAP2K1', 'agent=MAPK1',
                                type='Phosphorylation')

    def test_query_with_other(self):
        """Test that we can get an ActiveForm."""
        self.__check_good_query(agent='MAPK1', type='ActiveForm')

    def test_bad_camel(self):
        """Test that a type can be poorly formatted and resolve correctly."""
        self.__check_good_query(agent='MAPK1', type='acTivefOrm')

    def test_big_query(self):
        """Load-test with several big queries."""
        self.__check_good_query(agent='AKT1', check_stmts=False)
        self.__check_good_query(agent='MAPK1', check_stmts=False)
        self.__check_good_query(agent='TP53', check_stmts=False)

    def test_query_with_hgnc_ns(self):
        """Test specifying HGNC as a namespace."""
        self.__check_good_query(subject='6871@HGNC', object='MAP2K1',
                                type='Phosphorylation')

    def test_query_with_text_ns(self):
        """Test specifying TEXT as a namespace."""
        self.__check_good_query(subject='ERK@TEXT', type='Phosphorylation')

    def test_query_with_hgnc_symbol_ns(self):
        """Test specifying HGNC-SYMBOL as a namespace."""
        self.__check_good_query(subject='MAPK1@HGNC-SYMBOL',
                                type='Phosphorylation')

    def test_query_with_chebi_ns(self):
        """Test specifying CHEBI as a namespace."""
        self.__check_good_query(subject='CHEBI:6801@CHEBI')

    def test_query_with_bad_hgnc(self):
        resp, dt, size = self.__time_get_query('subject=MEK&object=ERK'
                                               '&type=Phosphorylation')
        assert resp.status_code != 200, "Got good status code."
        assert dt <= TIMELIMIT, dt
        assert size <= SIZELIMIT, size


if __name__ == '__main__':
    unittest.main()

