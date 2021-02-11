class AdaptiveAssembler:
    def __init__(self, unique_statements, filters, matches_fun=None):
        self.filters = filters
        self.matches_fun = matches_fun
        self.unique_statements = unique_statements
        self.stmts_by_hash = {stmt.get_hash(matches_fun=matches_fun): stmt
                              for stmt in self.unique_statements}
        assert len(self.stmts_by_hash) == len(self.unique_statements)
        for filter in self.filters:
            filter.initialize(self.stmts_by_hash)

    def get_all_refinements(self):
        all_refinements = []
        for sh, stmt in self.stmts_by_hash.items():
            all_refinements += [(sh, ref) for ref in
                                self.get_less_specifics(stmt)]
        return all_refinements

    def get_more_specifics(self, stmt):
        possibly_related = None
        for filter in self.filters:
            possibly_related = \
                filter.get_more_specifics(
                    stmt, possibly_related=possibly_related)
        return possibly_related

    def get_less_specifics(self, stmt):
        possibly_related = None
        for filter in self.filters:
            possibly_related = \
                filter.get_less_specifics(
                    stmt, possibly_related=possibly_related)
        return possibly_related
