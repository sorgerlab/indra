from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import numpy
import networkx
import indra.preassembler.sitemapper as sm


class BeliefEngine(object):
    """Assigns beliefs to INDRA Statements based on supporting evidence."""
    def __init__(self):
        self.prior_probs = {
            'rand': {
                'biopax': 0.2,
                'bel': 0.1,
                'trips': 0.3,
                'reach': 0.3,
                'biogrid': 0.01,
                'sparser': 0.3,
                'r3': 0.1,
                'phosphosite': 0.01,
                'assertion': 0.0
                },
            'syst': {
                'biopax': 0.01,
                'bel': 0.01,
                'trips': 0.05,
                'reach': 0.05,
                'biogrid': 0.01,
                'sparser': 0.3,
                'r3': 0.1,
                'phosphosite': 0.01,
                'assertion': 0.0
                }
            }

    def set_prior_probs(self, statements):
        """Sets the prior belief probabilities for a list of INDRA Statements.

        The Statements are assumed to be de-duplicated. In other words,
        each Statement in the list passed to this function is assumed to have
        a list of Evidence objects that support it. The prior probability of
        each Statement is calculated based on the number of Evidences it has
        and their sources.

        Parameters
        ----------
        statements : list[indra.statements.Statement]
            A list of INDRA Statements whose belief scores are to
            be calculated. Each Statement object's belief attribute is updated
            by this function.
        """
        for st in statements:
            sources = [ev.source_api for ev in st.evidence]
            uniq_sources = numpy.unique(sources)
            syst_factors = {s: self.prior_probs['syst'][s]
                            for s in uniq_sources}
            rand_factors = {k: [] for k in uniq_sources}
            for s in sources:
                rand_factors[s].append(self.prior_probs['rand'][s])
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
        """Sets hierarchical belief probabilities for a list of INDRA Statements.

        The Statements are assumed to be in a hierarchical relation graph with
        the supports and supported_by attribute of each Statement object having
        been set.
        The hierarchical belief probability of each Statement is calculated
        based on its prior probability and the probabilities propagated from
        Statements supporting it in the hierarchy graph.

        Parameters
        ----------
        statements : list[indra.statements.Statement]
            A list of INDRA Statements whose belief scores are to
            be calculated. Each Statement object's belief attribute is updated
            by this function.
        """
        def build_hierarchy_graph(stmts):
            g = networkx.DiGraph()
            for st1 in stmts:
                g.add_node(st1.matches_key(), {'stmt': st1})
                for st2 in st1.supported_by:
                    g.add_node(st2.matches_key(), {'stmt': st2})
                    g.add_edge(st2.matches_key(), st1.matches_key())
            return g
        def get_ranked_stmts(g):
            node_ranks = networkx.topological_sort(g, reverse=True)
            stmts = [g.node[n]['stmt'] for n in node_ranks]
            return stmts
        g = build_hierarchy_graph(statements)
        ranked_stmts = get_ranked_stmts(g)
        new_beliefs = []
        for st in ranked_stmts:
            bps = _get_belief_package(st)
            beliefs = [bp[0] for bp in bps]
            belief = 1 - numpy.prod([(1-b) for b in beliefs])
            new_beliefs.append(belief)
        for st, bel in zip(ranked_stmts, new_beliefs):
            st.belief = bel

    def set_linked_probs(self, linked_statements):
        """Sets the belief probabilities for a list of linked INDRA Statements.

        The list of LinkedStatement objects is assumed to come from the
        MechanismLinker. The belief probability of the inferred Statement is
        assigned the joint probability of its source Statements.

        Parameters
        ----------
        linked_statements : list[indra.mechlinker.LinkedStatement]
            A list of INDRA LinkedStatements whose belief scores are to
            be calculated. The belief attribute of the inferred Statement in
            the LinkedStatement object is updated by this function.
        """
        for st in linked_statements:
            source_probs = [s.belief for s in st.source_stmts]
            st.inferred_stmt.belief = numpy.prod(source_probs)


def _get_belief_package(stmt, n=1):
    def belief_stmts(belief_pkgs):
        return [pkg[1] for pkg in belief_pkgs]

    belief_packages = []
    for st in stmt.supports:
        parent_packages = _get_belief_package(st, n+1)
        belief_st = belief_stmts(belief_packages)
        for package in parent_packages:
            if not package[1] in belief_st:
                belief_packages.append(package)

    belief_package = (stmt.belief, stmt.matches_key())
    belief_packages.append(belief_package)
    return belief_packages

