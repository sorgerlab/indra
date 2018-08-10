import unittest
import json
import sys

from itertools import combinations
from datetime import datetime
from warnings import warn

from db_rest_api.api import MAX_STATEMENTS
from indra.db import get_primary_db
from indra.statements import stmts_from_json
from indra.databases import hgnc_client

from db_rest_api import api


TIMEGOAL = 1
TIMELIMIT = 30
SIZELIMIT = 4e7


def _check_stmt_agents(resp, agents):
    json_stmts = json.loads(resp.data.decode('utf-8'))['statements']
    stmts = stmts_from_json(json_stmts)
    for stmt in stmts:
        for ag_ix, db_ns, db_id in agents:
            if ag_ix is not None:
                assert stmt.agent_list()[ag_ix].db_refs.get(db_ns) == db_id
            # If the ag_ix is None, we're just checking that the Statement
            # contains the agent we're looking for
            else:
                db_ids = [ag.db_refs.get(db_ns) for ag in stmt.agent_list()]
                assert db_id in db_ids


class TimeWarning(Warning):
    pass


class DbApiTestCase(unittest.TestCase):

    def setUp(self):
        api.app.testing = True
        self.app = api.app.test_client()

    def tearDown(self):
        pass

    def __check_time(self, dt, time_goal=TIMEGOAL):
        print(dt)
        assert dt <= TIMELIMIT, \
            ("Query took %f seconds. Must be less than %f seconds."
             % (dt, TIMELIMIT))
        if dt >= TIMEGOAL:
            warn("Query took %f seconds, goal is less than %f seconds."
                 % (dt, TIMEGOAL), TimeWarning)
        return

    def __time_get_query(self, end_point, query_str):
        return self.__time_query('get', end_point, query_str)

    def __time_query(self, method, end_point, query_str=None, **data):
        start_time = datetime.now()
        if query_str is not None:
            url = '/%s/?%s' % (end_point, query_str)
        else:
            url = end_point
        meth_func = getattr(self.app, method)
        if data:
            resp = meth_func(url, data=json.dumps(data),
                             headers={'content-type': 'application/json'})
        else:
            resp = meth_func(url)
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
        time_goal = max(kwargs.pop('time_goal', TIMEGOAL), TIMEGOAL)
        query_str = '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()]
                             + list(args))
        resp, dt, size = self.__time_get_query('statements', query_str)
        assert resp.status_code == 200, \
            ('Got error code %d: \"%s\".'
             % (resp.status_code, resp.data.decode()))
        resp_dict = json.loads(resp.data.decode('utf-8'))
        #assert not resp_dict['limited']
        assert size <= SIZELIMIT, \
            ("Query took up %f MB. Must be less than %f MB."
             % (size/1e6, SIZELIMIT/1e6))
        self.__check_stmts(resp_dict['statements'], check_stmts=check_stmts)

        self.__check_time(dt, time_goal)
        return resp

    def __check_stmts(self, json_stmts, check_support=False, check_stmts=False):
        assert len(json_stmts) is not 0, \
            'Did not get any statements.'
        stmts = stmts_from_json(json_stmts)
        assert all([s.evidence for s in stmts]), \
            "Some statements lack evidence."

        # To allow for faster response-times, we currently do not include
        # support links in the response.
        if check_support:
            assert any([s.supports + s.supported_by for s in stmts]),\
                ("Some statements lack support: %s."
                 % str([str(s) for s in stmts if not s.supports+s.supported_by]))
            if check_stmts:
                assert all([not s1.matches(s2)
                            for s1, s2 in combinations(stmts, 2)]),\
                    ("Some statements match: %s."
                     % str([(s1, s2) for s1, s2 in combinations(stmts, 2)
                            if s1.matches(s2)]))
        return

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
        resp = self.__check_good_statement_query(object='MAP2K1',
                                                 subject='MAPK1',
                                                 type='Phosphorylation')
        _check_stmt_agents(resp, agents=[
                (0, 'HGNC', hgnc_client.get_hgnc_id('MAPK1')),
                (1, 'HGNC', hgnc_client.get_hgnc_id('MAP2K1'))])

    def test_object_only_query(self):
        """Test whether we can get an object only statement."""
        resp = self.__check_good_statement_query(object='GLUL',
                                          type='IncreaseAmount')
        _check_stmt_agents(resp, agents=[
                (1, 'HGNC', hgnc_client.get_hgnc_id('GLUL'))])
        return

    def test_query_with_two_agents(self):
        """Test a query were the roles of the agents are not given."""
        resp = self.__check_good_statement_query('agent=MAP2K1', 'agent=MAPK1',
                                                 type='Phosphorylation')
        _check_stmt_agents(resp, agents=[
                (None, 'HGNC', hgnc_client.get_hgnc_id('MAPK1')),
                (None, 'HGNC', hgnc_client.get_hgnc_id('MAP2K1'))])
        return

    def test_query_with_other(self):
        """Test that we can get an ActiveForm."""
        resp = self.__check_good_statement_query(agent='MAPK1',
                                                 type='ActiveForm')
        _check_stmt_agents(resp, agents=[
                (0, 'HGNC', hgnc_client.get_hgnc_id('MAPK1'))])
        return

    def test_bad_camel(self):
        """Test that a type can be poorly formatted and resolve correctly."""
        resp = self.__check_good_statement_query(agent='MAPK1',
                                                 type='acTivefOrm')
        _check_stmt_agents(resp, agents=[
                (0, 'HGNC', hgnc_client.get_hgnc_id('MAPK1'))])
        return

    def test_big_query(self):
        """Load-test with several big queries."""
        # Note that in this test we do not check the quality of the statements,
        # because there are likely to be so many statements that that would take
        # longer than needed, given that the quality is tested in other tests.
        self.__check_good_statement_query(agent='AKT1', check_stmts=False,
                                          time_goal=10)
        self.__check_good_statement_query(agent='MAPK1', check_stmts=False,
                                          time_goal=20)
        return

    def test_query_with_too_many_stmts(self):
        """Test our check of statement length and the response."""
        # As with the above, we do not check the quality of statements, because
        # there will be in general so many that it would make this test take far
        # to long, especially given that the statement quality is not really in
        # the scope of this test.
        resp, dt, size = self.__time_get_query('statements',
                                               'agent=TP53&on_limit=error')
        # assert resp.status_code == 413, "Unexpected status code: %s" % str(resp)
        assert dt < TIMELIMIT, "Query took too long: %s" % dt
        #assert 'Acetylation' in json.loads(resp.data.decode('utf-8'))\
        #                                  ['statements']
        resp, dt, size = self.__time_get_query('statements',
                                               'agent=TP53&on_limit=sample')
        # assert resp.status_code == 200, str(resp)
        assert dt < TIMELIMIT, dt
        resp_dict = json.loads(resp.data.decode('utf-8'))
        #assert len(resp_dict['statements']) == MAX_STATEMENTS
        resp, dt, size = self.__time_get_query('statements',
                                               'agent=TP53&on_limit=truncate')
        assert dt < TIMELIMIT, dt
        resp_dict = json.loads(resp.data.decode('utf-8'))
        #assert len(resp_dict['statements']) == MAX_STATEMENTS
        return

    def test_query_with_hgnc_ns(self):
        """Test specifying HGNC as a namespace."""
        resp = self.__check_good_statement_query(subject='6871@HGNC',
                                                 object='MAP2K1',
                                                 type='Phosphorylation')
        _check_stmt_agents(resp, agents=[
                (0, 'HGNC', '6871'),
                (1, 'HGNC', hgnc_client.get_hgnc_id('MAP2K1'))])
        return

    def test_query_with_text_ns(self):
        """Test specifying TEXT as a namespace."""
        resp = self.__check_good_statement_query(subject='ERK@TEXT',
                                          type='Phosphorylation')
        _check_stmt_agents(resp, agents=[(0, 'TEXT', 'ERK')])
        return

    def test_query_with_hgnc_symbol_ns(self):
        """Test specifying HGNC-SYMBOL as a namespace."""
        resp = self.__check_good_statement_query(subject='MAPK1@HGNC-SYMBOL',
                                                 type='Phosphorylation')
        _check_stmt_agents(resp, agents=[
                (0, 'HGNC', hgnc_client.get_hgnc_id('MAPK1'))])
        return

    def test_query_with_chebi_ns(self):
        """Test specifying CHEBI as a namespace."""
        resp = self.__check_good_statement_query(subject='CHEBI:6801@CHEBI')
        _check_stmt_agents(resp, agents=[(0, 'CHEBI', 'CHEBI:6801')])
        return

    def test_query_with_bad_hgnc(self):
        resp, dt, size = self.__time_get_query('statements',
                                               ('subject=MEK&object=ERK'
                                                '&type=Phosphorylation'))
        assert resp.status_code != 200, "Got good status code."
        self.__check_time(dt)
        assert size <= SIZELIMIT, size

    def test_famplex_query(self):
        resp, dt, size = self.__time_get_query('statements',
                                               ('object=PPP1C@FPLX'
                                                '&subject=CHEBI:44658@CHEBI'
                                                '&type=Inhibition'))
        resp_dict = json.loads(resp.data.decode('utf-8'))
        stmts = stmts_from_json(resp_dict['statements'])
        assert len(stmts)
        _check_stmt_agents(resp, agents=[
                (0, 'CHEBI', 'CHEBI:4465'),
                (1, 'FPLX', 'PPP1C')])
        self.__check_time(dt)
        assert size <= SIZELIMIT, size

    def test_statements_by_hashes_query(self):
        resp, dt, size = self.__time_query('get', 'statements/from_hashes',
                                           hashes=[-36028793042562873,
                                                   -12978096432588272,
                                                   -12724735151233845])
        resp_dict = json.loads(resp.data.decode('utf-8'))
        self.__check_stmts(resp_dict['statements'], check_support=True)
        self.__check_time(dt)
        return

    def test_statements_by_hashes_large_query(self):
        # Get a set of hashes.
        db = get_primary_db()
        res = db.select_sample_from_table(1000, db.EvidenceCounts)
        hash_cnt_dict = {ev_cts.mk_hash: ev_cts.ev_count for ev_cts in res}

        # Run the test.
        resp, dt, size = self.__time_query('get', 'statements/from_hashes',
                                           hashes=list(hash_cnt_dict.keys()))
        resp_dict = json.loads(resp.data.decode('utf-8'))
        self.__check_stmts(resp_dict['statements'], check_support=True)
        self.__check_time(dt, time_goal=20)
        return

    def __test_basic_paper_query(self, id_val, id_type, min_num_results=1):
        query_str = 'id=%s&type=%s' % (id_val, id_type)
        resp, dt, size = self.__time_get_query('papers', query_str)
        self.__check_time(dt)
        assert size <= SIZELIMIT, size
        assert resp.status_code == 200, str(resp)
        json_str = resp.data.decode('utf-8')
        json_list = json.loads(json_str)['statements']
        assert len(json_list) >= min_num_results, (min_num_results,
                                                   len(json_list))
        return

    def test_pmid_paper_query(self):
        pmid = '27014235'
        self.__test_basic_paper_query(pmid, 'pmid')

        # Now check without pmid specified (should be assumed.)
        resp, _, _ = self.__time_get_query('papers', 'id=%s' % pmid)
        assert resp.status_code == 200, str(resp)

    def test_pmcid_paper_query(self):
        self.__test_basic_paper_query('PMC5770457', 'pmcid')

    def test_trid_paper_query(self):
        self.__test_basic_paper_query('19649148', 'trid')


if __name__ == '__main__':
    unittest.main()

