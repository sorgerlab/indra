import unittest
from datetime import datetime
from time import sleep
from unittest import SkipTest

import pytest
from indra.sources import indra_db_rest as dbr
from indra.sources.indra_db_rest.api import get_statement_queries
from indra.sources.indra_db_rest.query import HasAgent, HasEvidenceBound
from indra.statements import Agent, Phosphorylation


EXPECTED_BATCH_SIZE = 500


def __check_request(seconds, *args, **kwargs):
    check_stmts = kwargs.pop('check_stmts', True)
    now = datetime.now()
    resp = dbr.get_statements(*args, **kwargs)
    time_taken = datetime.now() - now
    if check_stmts:
        assert resp.statements, "Got no statements."
    return resp


@pytest.mark.nonpublic
def test_simple_request():
    __check_request(6, 'MAP2K1', 'MAPK1', stmt_type='Phosphorylation')


@pytest.mark.nonpublic
def test_request_for_complex():
    __check_request(30, agents=['MEK@FPLX', 'ERK@FPLX'], stmt_type='Complex')


@pytest.mark.nonpublic
def test_null_request():
    try:
        dbr.get_statements()
    except ValueError:
        return
    except BaseException as e:
        assert False, "Raised wrong exception: " + str(e)
    assert False, "Null request did not raise any exception."


@pytest.mark.nonpublic
@pytest.mark.slow
def test_large_request():
    __check_request(40, agents=['AKT1'])


@pytest.mark.nonpublic
@pytest.mark.slow
def test_bigger_request():
    __check_request(60, agents=['MAPK1'])


@pytest.mark.nonpublic
def test_timeout_no_persist_gcg():
    resp = dbr.get_statements(agents=["GCG"],
                              persist=False,
                              timeout=0,)
    assert resp.is_working(), "Lookup resolved too fast."
    resp.wait_until_done(40)  # typically 20-30 s when slow/uncached
    assert len(resp.statements) > 0.9*EXPECTED_BATCH_SIZE, len(resp.statements)


@pytest.mark.nonpublic
def test_timeout_no_persist_type_object():
    resp = dbr.get_statements(stmt_type='phosphorylation',
                              object="NFkappaB@FPLX",
                              persist=False,
                              timeout=0)
    assert resp.is_working(), "Lookup resolved too fast."
    resp.wait_until_done(70)
    assert len(resp.statements) > 0.9*EXPECTED_BATCH_SIZE, len(resp.statements)


@pytest.mark.nonpublic
@pytest.mark.slow
def test_too_big_request_no_persist():
    resp_some = __check_request(60, agents=['TP53'], persist=False)
    assert sum(resp_some.get_ev_count(s) is not None
               for s in resp_some.statements) == len(resp_some.statements), \
        "Counts dict was improperly handled."
    return resp_some


@pytest.mark.nonpublic
@pytest.mark.slow
@pytest.mark.nogha
@unittest.skip('skipping')
def test_too_big_request_persist_and_block():
    resp_all1 = __check_request(200, agents=['TP53'], persist=True,
                                timeout=None)
    assert sum(resp_all1.get_ev_count(s) is not None
               for s in resp_all1.statements) > 0.9*len(resp_all1.statements), \
        "Counts dict was improperly handled."
    return resp_all1


@pytest.mark.nonpublic
@pytest.mark.slow
@pytest.mark.nogha
def test_too_big_request_persist_no_block():
    resp_some = test_too_big_request_no_persist()
    resp_all1 = test_too_big_request_persist_and_block()
    resp_all2 = __check_request(60, agents=['TP53'], persist=True,
                                timeout=10, check_stmts=False)
    num_counts = sum(resp_all2.get_ev_count(s) is not None
                     for s in resp_all2.statements)
    num_stmts = len(resp_all2.statements)
    assert num_counts == num_stmts, \
        ("Counts dict was improperly handled before completing: %d counts "
         "for %d statements." % (num_counts, num_stmts))
    assert resp_all2.is_working(), "Background complete resolved too fast."
    assert len(resp_all2.statements_sample) == len(resp_some.statements), \
        "Sample size: %s, Small resp size: %s" \
        % (len(resp_all2.statements_sample), len(resp_some.statements))
    resp_all2.wait_until_done(500)
    assert not resp_all2.is_working(), \
        "Response is still working. Took too long."
    assert resp_all2.statements
    num_counts = sum(resp_all2.get_ev_count(s) is not None
                     for s in resp_all2.statements)
    num_stmts = len(resp_all2.statements)
    assert num_counts > 0.9*num_stmts, \
        "Counts dict was improperly handled after completing."
    assert len(resp_all2.statements) == len(resp_all1.statements), \
        'Expected: %d, actual: %d' % (len(resp_all1.statements),
                                      len(resp_all2.statements))
    return


@pytest.mark.nonpublic
def test_famplex_namespace():
    p = dbr.get_statements('PDGF@FPLX', 'FOS', stmt_type='IncreaseAmount')
    stmts = p.statements
    print(len(stmts))
    print(stmts)
    assert all([s.agent_list()[0].db_refs.get('FPLX') == 'PDGF' for s in stmts]),\
        'Not all subjects match.'
    assert all([s.agent_list()[1].name == 'FOS' for s in stmts]),\
        'Not all objects match: ' \
        + ', '.join({s.agent_list()[1].name for s in stmts})


@pytest.mark.nonpublic
@pytest.mark.nogha
def test_paper_query():
    p = dbr.get_statements_for_papers([('pmcid', 'PMC5770457'),
                                       ('pmid', '27014235')])
    stmts_1 = p.statements
    assert len(stmts_1)

    p = dbr.get_statements_for_papers([('pmcid', 'PMC5770457'),
                                       ('pmid', '27014235')])
    assert len(p.statements)
    assert len(p.get_source_counts())
    assert len(p.get_ev_counts())


@pytest.mark.nonpublic
def test_regulate_amount():
    idbp = dbr.get_statements('FOS', stmt_type='RegulateAmount')
    stmts = idbp.statements
    stmt_types = {type(s).__name__ for s in stmts}
    counts = idbp.get_source_counts()
    one_key = list(counts.keys())[0]
    assert counts[one_key] == idbp.get_source_count_by_hash(one_key)
    assert {'IncreaseAmount', 'DecreaseAmount'}.issubset(stmt_types), \
        stmt_types


@pytest.mark.nonpublic
def test_get_statements_by_hash():
    hash_list = [30674674032092136, -22289282229858243, -25056605420392180]
    p = dbr.get_statements_by_hash(hash_list)
    stmts = p.statements
    print({s.get_hash(shallow=True): s for s in stmts})
    assert len(stmts) >= 2, \
        "Wrong number of statements: %s vs. %s" % (len(stmts), len(hash_list))

    p = dbr.get_statements_by_hash(hash_list)
    assert len(p.statements)
    assert len(p.get_source_counts())
    assert len(p.get_ev_counts())
    return


@pytest.mark.nonpublic
def test_get_statements_by_hash_no_hash():
    p = dbr.get_statements_by_hash([])
    assert not p.statements, "Got statements without giving a hash."


@pytest.mark.nonpublic
def test_curation_submission():
    from indra.config import get_config
    api_key = get_config('INDRA_DB_REST_API_KEY', failure_ok=True)
    if not api_key:
        raise SkipTest("No API Key, this test will not work.")
    res = dbr.submit_curation(32760831642168299, 'TEST', 'This is a test.',
                              'tester', is_test=True)
    assert res['result'] == 'test passed', res


@pytest.mark.nonpublic
def test_get_curations():
    from indra.config import get_config
    api_key = get_config('INDRA_DB_REST_API_KEY', failure_ok=True)
    if not api_key:
        raise SkipTest("No API Key, this test will not work.")
    stmt_hash = -13159234982749425
    src_hash = 3817298406742073624
    res = dbr.get_curations(stmt_hash, src_hash)
    assert isinstance(res, list)
    assert len(res) > 0
    assert all(c['pa_hash'] == stmt_hash and c['source_hash'] == src_hash
               for c in res)


@pytest.mark.nonpublic
def test_get_statement_queries():
    ag = Agent('MAP2K1', db_refs={})
    stmt = Phosphorylation(None, ag)
    urls = get_statement_queries([stmt])
    assert 'MAP2K1@NAME' in urls[0]
    urls = get_statement_queries([stmt], fallback_ns='TEXT')
    assert 'MAP2K1@TEXT' in urls[0]
    urls = get_statement_queries([stmt],
                                 pick_ns_fun=lambda x: '%s@%s' %
                                                       (x.name, 'XXX'))
    assert 'MAP2K1@XXX' in urls[0], urls[0]
    ag = Agent('MEK', db_refs={'FPLX': 'MEK'})
    stmt = Phosphorylation(None, ag)
    urls = get_statement_queries([stmt])
    assert 'MEK@FPLX' in urls[0]
    urls = get_statement_queries([stmt], fallback_ns='TEXT')
    assert 'MEK@FPLX' in urls[0]
    urls = get_statement_queries([stmt],
                                 pick_ns_fun=lambda x: '%s@%s' %
                                                       (x.name, 'XXX'))


@pytest.mark.nonpublic
def test_get_statements_end_on_limit():
    p = dbr.get_statements(subject="TNF", limit=1400, timeout=1)
    try:
        t = 0
        violations = 0
        violations_allowed = 3
        while p.is_working():
            assert t < 100, t
            limit = p._get_next_limit()
            if limit == 0 and p.is_working():
                violations += 1
                assert violations <= violations_allowed
            sleep(1)
            t += 1
    finally:
        p.cancel()
        p.wait_until_done()


@pytest.mark.nonpublic
def test_get_statements_evidence_bounded():
    query = HasAgent('MEK') & HasEvidenceBound(["< 10"])
    p = dbr.get_statements_from_query(query, limit=10)
    stmts = p.statements
    assert len(stmts) == 10
    assert all(c < 10 for c in p.get_ev_counts().values())


@pytest.mark.nonpublic
def test_get_statements_strict_stop_short():
    start = datetime.now()
    p = dbr.get_statements("TNF", timeout=1, strict_stop=True)
    end = datetime.now()
    sleep(5)
    assert not p.is_working()
    dt = (end - start).total_seconds()
    assert 1 <= dt < 1.5, dt
    assert not p.statements
    assert not p.statements_sample


@pytest.mark.nonpublic
def test_get_statements_strict_stop_long():
    timeout = 40  # Typically 20-30 s when slow/uncached
    start = datetime.now()
    p = dbr.get_statements("VEGF", timeout=timeout, strict_stop=True)
    end = datetime.now()
    sleep(5)
    assert not p.is_working()
    dt = (end - start).total_seconds()
    assert timeout <= dt < (timeout + 0.5), dt
    assert p.statements


@pytest.mark.nonpublic
@pytest.mark.nogha
def test_filter_ev():
    ids = [('pmcid', 'PMC5770457'), ('pmid', '27014235')]
    p = dbr.get_statements_for_papers(ids)
    assert p.statements

    correct_source = 0
    incorrect_source = 0
    for s in p.statements:
        for ev in s.evidence:
            if any(ev.text_refs.get(t.upper()) == v for t, v in ids):
                correct_source += 1
            else:
                incorrect_source += 1

    assert incorrect_source == 0,\
        f"{incorrect_source} unfiltered sources vs. {correct_source} filtered."


@pytest.mark.nonpublic
def test_sort_by_belief():
    p = dbr.get_statements(object="MEK", stmt_type="Inhibition",
                                  sort_by='belief', limit=10)
    assert p.statements
    beliefs = [s.belief for s in p.statements]
    assert beliefs == sorted(beliefs, reverse=True), \
        f"Beliefs mis-ordered!\nbeliefs: {beliefs}\n" \
        f"belief_dict: {p.get_belief_scores()}"


@pytest.mark.nonpublic
def test_sort_by_ev_count():
    p = dbr.get_statements(object="MEK", stmt_type="Inhibition",
                           sort_by='ev_count', limit=10, ev_limit=None)
    assert p.statements
    counts = [p.get_ev_count(s) for s in p.statements]
    assert counts == sorted(counts, reverse=True),\
        f"Counts mis-ordered!\ncounts: {counts}\nev_counts: {p.get_ev_counts()}"


@pytest.mark.nonpublic
@unittest.skip('This query test fails intermittently')
def test_namespace_only_agent_query():
    q = HasAgent("MEK") & HasAgent(namespace="CHEBI")
    p = dbr.get_statements_from_query(q, limit=10)
    assert p.statements
    assert all(any("CHEBI" in ag.db_refs for ag in s.agent_list())
               and any(ag.db_refs.get("FPLX") == "MEK" for ag in s.agent_list())
               for s in p.statements)
