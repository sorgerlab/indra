from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from datetime import datetime

from nose.plugins.attrib import attr
from indra.sources import indra_db_rest as dbr
from indra.sources.indra_db_rest import IndraDBRestError


def __check_request(seconds, *args, **kwargs):
    check_stmts = kwargs.pop('check_stmts', True)
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
    __check_request(5, 'MAP2K1', 'MAPK1', stmt_type='Phosphorylation')


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
    __check_request(20, agents=['AKT1'])


@attr('nonpublic')
def test_bigger_request():
    __check_request(30, agents=['MAPK1'])


@attr('nonpublic')
def test_too_big_request():
    resp_some = __check_request(30, agents=['TP53'], persist=False,
                                simple_response=False)
    resp_all1 = __check_request(60, agents=['TP53'], persist=True, block=True,
                                simple_response=False)
    resp_all2 = __check_request(30, agents=['TP53'], persist=True, block=False,
                                check_stmts=False, simple_response=False)
    assert not resp_all2.done, \
        "Background complete appears to have resolved too fast."
    assert len(resp_all2.statements_sample) == len(resp_some.statements)
    resp_all2.wait_until_done(60)
    assert len(resp_all2.statements) == len(resp_all1.statements)
    return


@attr('nonpublic')
def test_famplex_namespace():
    stmts = dbr.get_statements('PDGF@FPLX', 'FOS', stmt_type='IncreaseAmount')
    print(len(stmts))
    assert all([s.agent_list()[0].db_refs.get('FPLX') == 'PDGF' for s in stmts]),\
        'Not all subjects match.'
    assert all([s.agent_list()[1].name == 'FOS' for s in stmts]),\
        'Not all objects match.'


@attr('nonpublic')
def test_paper_query():
    stmts_1 = dbr.get_statements_for_paper('PMC495400', 'pmcid')
    assert len(stmts_1)
    stmts_2 = dbr.get_statements_for_paper('14461663')
    assert len(stmts_2)


@attr('nonpublic')
def test_regulate_amount():
    stmts = dbr.get_statements('FOS', stmt_type='RegulateAmount')
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
    try:
        stmts = dbr.get_statements_by_hash([])
    except IndraDBRestError as e:
        assert e.status_code == 400, \
            "Query failed for wrong reason:\n%s" % str(e)
        return
    assert False, "Query with no hashes did not get an exception."
