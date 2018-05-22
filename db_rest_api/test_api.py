import unittest
import json
import sys

from itertools import combinations
from datetime import datetime

from db_rest_api.api import MAX_STATEMENTS
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

    def __time_get_query(self, end_point, query_str):
        start_time = datetime.now()
        resp = self.app.get('/%s/?%s' % (end_point, query_str))
        t_delta = datetime.now() - start_time
        dt = t_delta.seconds + t_delta.microseconds/1e6
        print(dt)
        size = int(resp.headers['Content-Length'])
        raw_size = sys.getsizeof(resp.data)
        print("Raw size: {raw:f}/{lim:f}, Compressed size: {comp:f}/{lim:f}."
              .format(raw=raw_size/1e6, lim=SIZELIMIT/1e6, comp=size/1e6))
        return resp, dt, size

    def __check_good_statement_query(self, *args, **kwargs):
        check_stmts = kwargs.pop('check_stmts', True)
        time_limit = kwargs.pop('time_limit', TIMELIMIT)
        query_str = '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()]
                             + list(args))
        resp, dt, size = self.__time_get_query('statements', query_str)
        assert resp.status_code == 200, \
            ('Got error code %d: \"%s\".'
             % (resp.status_code, resp.data.decode()))
        resp_dict = json.loads(resp.data.decode('utf-8'))
        assert not resp_dict['limited']
        json_stmts = resp_dict['statements']
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
        resp, dt, size = self.__time_get_query('statements', '')
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
        self.__check_good_statement_query(object='MAP2K1', subject='MAPK1',
                                          type='Phosphorylation')

    def test_object_only_query(self):
        """Test whether we can get an object only statement."""
        self.__check_good_statement_query(object='GLUL',
                                          type='IncreaseAmount')

    def test_query_with_two_agents(self):
        """Test a query were the roles of the agents are not given."""
        self.__check_good_statement_query('agent=MAP2K1', 'agent=MAPK1',
                                          type='Phosphorylation')

    def test_query_with_other(self):
        """Test that we can get an ActiveForm."""
        self.__check_good_statement_query(agent='MAPK1', type='ActiveForm')

    def test_bad_camel(self):
        """Test that a type can be poorly formatted and resolve correctly."""
        self.__check_good_statement_query(agent='MAPK1', type='acTivefOrm')

    def test_big_query(self):
        """Load-test with several big queries."""
        self.__check_good_statement_query(agent='AKT1', check_stmts=False,
                                          time_limit=5)
        self.__check_good_statement_query(agent='MAPK1', check_stmts=False,
                                          time_limit=10)

    def test_query_with_too_many_stmts(self):
        """Test our check of statement length and the response."""
        resp, dt, size = self.__time_get_query('statements',
                                               'agent=TP53&on_limit=error')
        assert resp.status_code == 413, "Unexpected status code: %s" % str(resp)
        assert dt < 30, "Query took too long: %d" % dt
        assert 'Acetylation' in json.loads(resp.data.decode('utf-8'))['statements']
        resp, dt, size = self.__time_get_query('statements',
                                               'agent=TP53&on_limit=sample')
        assert resp.status_code == 200, str(resp)
        assert dt < 30, dt
        resp_dict = json.loads(resp.data.decode('utf-8'))
        assert len(resp_dict['statements']) == MAX_STATEMENTS
        resp, dt, size = self.__time_get_query('statements',
                                               'agent=TP53&on_limit=truncate')

    def test_query_with_hgnc_ns(self):
        """Test specifying HGNC as a namespace."""
        self.__check_good_statement_query(subject='6871@HGNC', object='MAP2K1',
                                          type='Phosphorylation')

    def test_query_with_text_ns(self):
        """Test specifying TEXT as a namespace."""
        self.__check_good_statement_query(subject='ERK@TEXT', type='Phosphorylation')

    def test_query_with_hgnc_symbol_ns(self):
        """Test specifying HGNC-SYMBOL as a namespace."""
        self.__check_good_statement_query(subject='MAPK1@HGNC-SYMBOL',
                                          type='Phosphorylation')

    def test_query_with_chebi_ns(self):
        """Test specifying CHEBI as a namespace."""
        self.__check_good_statement_query(subject='CHEBI:6801@CHEBI')

    def test_query_with_bad_hgnc(self):
        resp, dt, size = self.__time_get_query('statements',
                                               ('subject=MEK&object=ERK'
                                                '&type=Phosphorylation'))
        assert resp.status_code != 200, "Got good status code."
        assert dt <= TIMELIMIT, dt
        assert size <= SIZELIMIT, size

    def test_famplex_query(self):
        resp, dt, size = self.__time_get_query('statements',
                                               ('subject=PDGF@FPLX'
                                                '&object=FOS'
                                                '&type=Phosphorylation'))
        resp_dict = json.loads(resp.data.decode('utf-8'))
        stmts = stmts_from_json(resp_dict['statements'])
        assert len(stmts)
        assert all([s.agent_list()[0].db_refs.get('FPLX') == 'PDGF'
                    for s in stmts]),\
            'Not all subjects match.'
        assert dt <= TIMELIMIT, dt
        assert size <= SIZELIMIT, size

    def __test_basic_paper_query(self, id_val, id_type, min_num_results=1):
        query_str = 'id=%s&type=%s' % (id_val, id_type)
        resp, dt, size = self.__time_get_query('papers', query_str)
        assert dt <= TIMELIMIT, dt
        assert size <= SIZELIMIT, size
        assert resp.status_code == 200, str(resp)
        json_str = resp.data.decode('utf-8')
        json_list = json.loads(json_str)['statements']
        assert len(json_list) >= min_num_results, (min_num_results,
                                                   len(json_list))
        return

    def test_pmid_paper_query(self):
        self.__test_basic_paper_query('8436299', 'pmid')

        # Now check without pmid specified (should be assumed.)
        resp, _, _ = self.__time_get_query('papers', 'id=8436299')
        assert resp.status_code == 200, str(resp)

    def test_pmcid_paper_query(self):
        self.__test_basic_paper_query('PMC5770457', 'pmcid')

    def test_trid_paper_query(self):
        self.__test_basic_paper_query('28145129', 'trid')


if __name__ == '__main__':
    unittest.main()

