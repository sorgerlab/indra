from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import numpy
import indra.preassembler.sitemapper as sm

class BeliefEngine(object):
    def __init__(self, statements):
        self.statements = statements
        self.prior_probs_rand = {
            'biopax': 0.2,
            'bel': 0.1,
            'trips': 0.4,
            'reach': 0.3,
            'biogrid': 0.01,
            'assertion': 0.0
            }
        self.prior_probs_syst = {
            'biopax': 0.01,
            'bel': 0.01,
            'trips': 0.2,
            'reach': 0.,
            'biogrid': 0.01,
            'assertion': 0.0
            }

    def set_prior_probs(self, statements):
        for st in statements:
            sources = [ev.source_api for ev in st.evidence]
            uniq_sources = numpy.unique(sources)
            syst_factors = {s: self.prior_probs_syst[s] for s in uniq_sources}
            rand_factors = {k: [] for k in uniq_sources}
            for s in sources:
                rand_factors[s].append(self.prior_probs_rand[s])
            neg_prob_prior = 1
            for s in uniq_sources:
                neg_prob_prior *= (syst_factors[s] +
                                   numpy.prod(rand_factors[s]))
            prob_prior = 1 - neg_prob_prior
            vs, _ = sm.default_mapper.map_sites([st])
            if not vs:
                prob_prior *= 0.05
            st.belief = prob_prior

    def set_hierarchy_probs(self, statements):
        for st in statements:
            prob = self.get_rolling_prob(st)
            st.belief = prob

    def set_linked_probs(self, linked_statements):
        for st in linked_statements:
            source_probs = [s.belief for s in st.source_stmts]
            st.inferred_stmt.belief = numpy.prod(source_probs)

    @staticmethod
    def get_rolling_prob(stmt):
        neg_prob_self = (1-stmt.belief)
        neg_probs_rolling = 1
        for st in stmt.supports:
            neg_probs_rolling *= (1-self.get_rolling_prob(st))
        return (1-neg_prob_self*neg_probs_rolling)

