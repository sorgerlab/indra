import copy
from indra.statements import Event


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


def get_more_generic_agent(agent, ontology):
    generic_agent = copy.deepcopy(agent)
    db_ns, db_id = agent.get_grounding()
    if db_ns is not None:
        parents = ontology.get_parents(db_ns, db_id)


def generate_generics(stmt):
    if isinstance(stmt, Event):
        generalized_event = copy.deepcopy(stmt)
        generalized_concept = \
            generalize_concept_grounding(generalized_event.concept)
        generalized_event.concept = generalized_concept

def generalize_concept_grounding(concept):
    if 'WM' in concept.db_refs:
        wm_grounding = concept.db_refs['WM']
