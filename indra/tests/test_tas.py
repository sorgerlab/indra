from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.sources.tas.api import _load_data, process_csv


def test_load_data():
    data = _load_data()
    assert len(data) > 100, len(data)


def test_processor():
    tp = process_csv()
    assert tp
    assert tp.statements
    num_stmts = len(tp.statements)
    num_data = len(_load_data())
    assert num_stmts == num_data, \
        "Expected %d stmts, got %d." % (num_data, num_stmts)
