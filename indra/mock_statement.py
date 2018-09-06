class MockStatement(object):
    def __init__(evidence, mk_hash, supports):
        self.evidence = evidence
        self.__mk_hash = mk_hash
        self.supports = supports
        self.belief = None

    def matches_key():
        return self.__mk_hash


class MockEvidence(object):
    def __init__(source_api, selected_annotations):
        self.source_api = source_api
        # For 
        #    biopax: 'source_sub_id'
        #    reach: 'found_by'
        #    geneways: 'actiontype'
        self.annotations = selected_annotations
