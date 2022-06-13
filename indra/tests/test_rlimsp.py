import os
import unittest
from indra.sources import rlimsp

@unittest.skip('RLIMS-P webservice is down')
def test_simple_usage():
    rp = rlimsp.process_from_webservice('PMC3717945')
    stmts = rp.statements
    assert len(stmts) == 33, len(stmts)
    for s in stmts:
        assert len(s.evidence) == 1, "Wrong amount of evidence."
        ev = s.evidence[0]
        assert ev.annotations, "Missing annotations."
        assert 'agents' in ev.annotations.keys()
        assert 'trigger' in ev.annotations.keys()


@unittest.skip('RLIMS-P webservice is down')
def test_ungrounded_endpoint_with_pmids():
    pmid_list = ['16403219', '22258404', '16961925', '22096607']
    stmts = []
    for pmid in pmid_list:
        rp = rlimsp.process_from_webservice(pmid, id_type='pmid')
        assert len(rp.statements) > 10, len(rp.statements)
        stmts.extend(rp.statements)
    assert len(stmts) == 394, len(stmts)
    return


def test_tyrosine_grounding():
    here = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(here, 'rlimsp_site.json')
    rp = rlimsp.process_jsonl_file(fname)
    assert len(rp.statements) == 1
    stmt = rp.statements[0]
    assert stmt.residue == 'Y'
    assert stmt.position == '705'
