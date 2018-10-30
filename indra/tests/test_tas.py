from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.sources.tas.api import _load_data, process_csv


def test_load_data():
    data = _load_data()
    assert len(data) > 100, len(data)


def test_processor():
    tp = process_csv(affinity_class_limit=10)
    assert tp
    assert tp.statements
    num_stmts = len(tp.statements)
    # This is the total number of statements about human genes
    assert num_stmts == 51722, num_stmts
    assert all(len(s.evidence) == 1 for s in tp.statements), \
        "Some statements lack evidence, or have extra evidence."
