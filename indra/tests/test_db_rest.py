from __future__ import absolute_import, print_function, unicode_literals

import random
from builtins import dict, str

from datetime import datetime

from nose.plugins.attrib import attr
from indra.sources import indra_db_rest as dbr


EXPECTED_BATCH_SIZE = 1000


def __check_request(seconds, *args, **kwargs):
    check_stmts = kwargs.pop('check_stmts', True)
    simple_response = kwargs.pop('simple_response', True)
    kwargs['simple_response'] = simple_response
    now = datetime.now()
    resp = dbr.get_statements(*args, **kwargs)
    time_taken = datetime.now() - now
    if check_stmts:
        if kwargs.get('simple_response', True):
            stmts = resp
        else:
            stmts = resp.statements
        assert stmts, "Got no statements."
    assert time_taken.seconds < seconds, time_taken.seconds
    return resp


@attr('nonpublic')
def test_simple_request():
    __check_request(6, 'MAP2K1', 'MAPK1', stmt_type='Phosphorylation')


@attr('nonpublic')
def test_request_for_complex():
    __check_request(30, agents=['MEK@FPLX', 'ERK@FPLX'], stmt_type='Complex')


@attr('nonpublic')
def test_null_request():
    try:
        dbr.get_statements()
    except ValueError:
        return
    except BaseException as e:
        assert False, "Raised wrong exception: " + str(e)
    assert False, "Null request did not raise any exception."


@attr('nonpublic')
def test_large_request():
    __check_request(25, agents=['AKT1'])


@attr('nonpublic')
def test_bigger_request():
    __check_request(30, agents=['MAPK1'])


@attr('nonpublic')
def test_timeout_no_persist_agent():
    candidates = ['TP53', 'NFkappaB@FPLX', 'AKT@FPLX']
    agent = random.choice(candidates)
    print(agent)
    resp = dbr.get_statements(agents=[agent], persist=False, timeout=0)
    assert resp.is_working(), "Lookup resolved too fast."
    resp.wait_until_done(70)
    assert len(resp.statements) == EXPECTED_BATCH_SIZE, len(resp.statements)


@attr('nonpublic')
def test_timeout_no_persist_type_object():
    candidates = ['TP53', 'NFkappaB@FPLX', 'AKT@FPLX']
    agent = random.choice(candidates)
    print(agent)
    resp = dbr.get_statements(stmt_type='phosphorylation', object=agent,
                              persist=False, timeout=0)
    assert resp.is_working(), "Lookup resolved too fast."
    resp.wait_until_done(70)
    assert len(resp.statements) == EXPECTED_BATCH_SIZE, len(resp.statements)


@attr('nonpublic')
def test_too_big_request_no_persist():
    resp_some = __check_request(60, agents=['TP53'], persist=False,
                                simple_response=False)
    assert sum(resp_some.get_ev_count(s) is not None
               for s in resp_some.statements) == len(resp_some.statements), \
        "Counts dict was improperly handled."
    return resp_some


@attr('nonpublic', 'slow')
def test_too_big_request_persist_and_block():
    resp_all1 = __check_request(200, agents=['TP53'], persist=True, timeout=None,
                                simple_response=False)
    assert sum(resp_all1.get_ev_count(s) is not None
               for s in resp_all1.statements) > 0.9*len(resp_all1.statements), \
        "Counts dict was improperly handled."
    return resp_all1


@attr('nonpublic', 'slow')
def test_too_big_request_persist_no_block():
    resp_some = test_too_big_request_no_persist()
    resp_all1 = test_too_big_request_persist_and_block()
    resp_all2 = __check_request(60, agents=['TP53'], persist=True,
                                timeout=10, check_stmts=False,
                                simple_response=False)
    num_counts = sum(resp_all2.get_ev_count(s) is not None
                     for s in resp_all2.statements)
    num_stmts = len(resp_all2.statements)
    assert num_counts == num_stmts, \
        ("Counts dict was improperly handled before completing: %d counts "
         "for %d statements." % (num_counts, num_stmts))
    assert resp_all2.is_working(), "Background complete resolved too fast."
    assert len(resp_all2.statements_sample) == len(resp_some.statements)
    resp_all2.wait_until_done(120)
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


@attr('nonpublic')
def test_famplex_namespace():
    stmts = dbr.get_statements('PDGF@FPLX', 'FOS', stmt_type='IncreaseAmount',
                               simple_response=True)
    print(len(stmts))
    assert all([s.agent_list()[0].db_refs.get('FPLX') == 'PDGF' for s in stmts]),\
        'Not all subjects match.'
    assert all([s.agent_list()[1].name == 'FOS' for s in stmts]),\
        'Not all objects match.'


@attr('nonpublic')
def test_paper_query():
    stmts_1 = dbr.get_statements_for_paper([('pmcid', 'PMC5770457'),
                                            ('pmid', '27014235')])
    assert len(stmts_1)


@attr('nonpublic')
def test_regulate_amount():
    stmts = dbr.get_statements('FOS', stmt_type='RegulateAmount',
                               simple_response=True)
    print(len(stmts))
    stmt_types = {type(s).__name__ for s in stmts}
    print(stmt_types)
    assert {'IncreaseAmount', 'DecreaseAmount'}.issubset(stmt_types), stmt_types


@attr('nonpublic')
def test_get_statements_by_hash():
    hash_list = [-36028793042562873, -12978096432588272, -12724735151233845]
    stmts = dbr.get_statements_by_hash(hash_list)
    print({s.get_hash(shallow=True): s for s in stmts})
    assert len(stmts) == len(hash_list), \
        "Wrong number of statements: %s vs. %s" % (len(stmts), len(hash_list))
    return


@attr('nonpublic')
def test_get_statements_by_hash_no_hash():
    stmts = dbr.get_statements_by_hash([])
    assert not stmts, "Got statements without giving a hash."


@attr('nonpublic')
def test_curation_submission():
    dbr.submit_curation(-36028793042562873, 'TEST', 'This is a test.',
                        'tester', is_test=True)
