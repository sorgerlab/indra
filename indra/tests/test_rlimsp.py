from indra.sources import rlimsp


def test_simple_usage():
    rp = rlimsp.process_pmc('PMC3717945')
    stmts = rp.statements
    assert len(stmts) == 6, len(stmts)


def test_ungrounded_usage():
    rp = rlimsp.process_pmc('PMC3717945', with_grounding=False)
    assert len(rp.statements) == 33, len(rp.statements)
