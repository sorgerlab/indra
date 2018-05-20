from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from datetime import datetime

from nose.plugins.attrib import attr
from indra.sources import indra_db_rest as dbr


def __check_request(seconds, *args, **kwargs):
    now = datetime.now()
    stmts = dbr.get_statements(*args, **kwargs)
    assert stmts, "Got no statements."
    time_taken = datetime.now() - now
    assert time_taken.seconds < seconds, time_taken.seconds


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
    __check_request(60, agents=['TP53'])


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
    stmts_1 = dbr.get_statements_for_paper('PMC5770457', 'pmcid')
    assert len(stmts_1)
    stmts_2 = dbr.get_statements_for_paper('8436299')
    assert len(stmts_2)


@attr('nonpublic')
def test_regulate_amount():
    stmts = dbr.get_statements('FOS', stmt_type='RegulateAmount')
    print(len(stmts))
    stmt_types = {type(s).__name__ for s in stmts}
    print(stmt_types)
    assert {'IncreaseAmount', 'DecreaseAmount'}.issubset(stmt_types), stmt_types
