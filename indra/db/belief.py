class MockStatement(object):
    """A class to imitate real INDRA Statements for calculating belief."""
    def __init__(self, evidence, mk_hash, supports):
        if isinstance(evidence, list):
            self.evidence = evidence
        else:
            self.evidence = [evidence]
        self.__mk_hash = mk_hash
        self.supports = supports
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
