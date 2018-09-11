class MockStatement(object):
    """A class to imitate real INDRA Statements for calculating belief."""
    def __init__(self, mk_hash, evidence, supports=None, supported_by=None):
        if isinstance(evidence, list):
            self.evidence = evidence
        else:
            self.evidence = [evidence]
        self.__mk_hash = mk_hash
        if supports:
            self.supports = supports
        else:
            self.supports = []
        if supported_by:
            self.supported_by = supported_by
        else:
            self.supported_by = []
        self.belief = None

    def matches_key(self):
        return self.__mk_hash


class MockEvidence(object):
    """A class to imitate real INDRA Evidence for calculating belief."""
    def __init__(self, source_api):
        self.source_api = source_api

        # Some annotations are used in indra.belief.tag_evidence_subtype.
        # TODO: optionally implement necessary annotations.
        self.annotations = {}


def populate_support(stmts, links):
    """Populate the supports supported_by lists of statements given links.

    Parameters
    ----------
    stmts : list[MockStatement/Statement]
        A list of objects with supports and supported_by attributes which are
        lists or equivalent.
    links : list[tuple]
        A list of pairs of hashes or matches_keys, where the first supports the
        second.
    """
    stmt_dict = {s.matches_key(): s for s in stmts}
    for supped_idx, supping_idx in links:
        stmt_dict[supping_idx].supports.append(stmt_dict[supped_idx])
        stmt_dict[supped_idx].supported_by.append(stmt_dict[supping_idx])
    return
