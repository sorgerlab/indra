from indra.sources import rlimsp


def test_simple_usage():
    rp = rlimsp.process_pmc('PMC3717945')
    stmts = rp.statements
    assert len(stmts) == 6, len(stmts)
    for s in stmts:
        assert len(s.evidence) == 1, "Wrong amount of evidence."
        ev = s.evidence[0]
        assert ev.annotations, "Missing annotations."
        assert 'agents' in ev.annotations.keys()
        assert 'trigger' in ev.annotations.keys()


def test_ungrounded_usage():
    rp = rlimsp.process_pmc('PMC3717945', with_grounding=False)
    assert len(rp.statements) == 33, len(rp.statements)
