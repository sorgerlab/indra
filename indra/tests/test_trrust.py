import pytest
from indra.sources import trrust
from indra.statements import RegulateAmount


@pytest.mark.slow
@pytest.mark.webservice
def test_process_from_web():
    tp = trrust.process_from_web()
    assert len(tp.statements) > 6200
    for stmt in tp.statements:
        assert isinstance(stmt, RegulateAmount)
        assert len(stmt.evidence) == 1
        assert stmt.obj.db_refs.get('HGNC'), stmt.obj.db_refs
        assert stmt.subj.db_refs.get('HGNC'), stmt.subj.db_refs
        assert stmt.evidence[0].source_api == 'trrust'
        assert stmt.evidence[0].pmid is not None
