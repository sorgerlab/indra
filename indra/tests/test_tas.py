import pytest
from indra.sources.tas import process_from_web


@pytest.mark.slow
def test_processor():
    tp = process_from_web(affinity_class_limit=10)
    assert tp
    assert tp.statements
    num_stmts = len(tp.statements)
    # This is the total number of statements about human genes
    assert num_stmts == 1123724, num_stmts
    assert all(len(s.evidence) >= 1 for s in tp.statements), \
        'Some statements lack any evidence'
