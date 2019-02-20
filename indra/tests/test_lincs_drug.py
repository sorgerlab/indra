from __future__ import absolute_import, print_function, unicode_literals

import unittest
from nose.plugins.attrib import attr
from indra.databases.lincs_client import get_drug_target_data
from indra.sources.lincs_drug import process_from_web


@attr('webservice')
@unittest.skip('LINCS web service very unreliable.')
def test_process_from_web():
    lincs_p = process_from_web()
    assert lincs_p is not None
    assert lincs_p.statements
    data_len = len(get_drug_target_data())
    num_stmts = len(lincs_p.statements)
    # Note that due to an erroneous entry in the HMS LINCS protein table,
    # one Statement is not extracted from the table, hence the condition
    # below, which should be kept as long as the error persists
    assert num_stmts >= data_len - 1, \
        ("Did not convert all statements: expected %d, got %d."
         % (data_len, num_stmts))
    assert all(len(s.evidence) > 0 for s in lincs_p.statements),\
        "Some statements lack evidence."
