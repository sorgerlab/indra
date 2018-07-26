from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import numpy
import networkx
import logging

try:
    from indra.sources.reach.processor import determine_reach_subtype
    use_reach_subtypes = True
except ImportError:
    use_reach_subtypes = False

logger = logging.getLogger("belief")

default_probs = {
    'rand': {
        'biopax': 0.2,
        'bel': 0.1,
        'trips': 0.3,
        'reach': 0.3,
        'isi': 0.3,
        'eidos': 0.3,
        'hume': 0.3,
        'cwms': 0.3,
        'sofia': 0.3,
        'biogrid': 0.01,
        'sparser': 0.3,
        'r3': 0.1,
        'phosphosite': 0.01,
        'ndex': 0.049,
        'signor': 0.049,
        'assertion': 0.0,
        },
    'syst': {
        'biopax': 0.01,
        'bel': 0.01,
        'trips': 0.05,
        'reach': 0.05,
        'isi': 0.05,
        'eidos': 0.05,
        'hume': 0.05,
        'cwms': 0.05,
        'sofia': 0.05,
        'biogrid': 0.01,
        'sparser': 0.05,
        'r3': 0.05,
        'phosphosite': 0.01,
        'ndex': 0.01,
        'signor': 0.01,
        'assertion': 0.0,
        }
    }


class BeliefEngine(object):
    """Assigns beliefs to INDRA Statements based on supporting evidence.

    Parameters
    ----------
    prior_probs : Optional[dict[dict]]
        A dictionary of prior probabilities used to override/extend the
        default ones. There are two types of prior probabilities: rand and syst
        corresponding to random error and systematic error rate for each
        knowledge source. The prior_probs dictionary has the general structure
        {'rand': {'s1': pr1, ..., 'sn': prn},
        'syst': {'s1': ps1, ..., 'sn': psn}}
        where 's1' ... 'sn' are names of input sources and pr1 ... prn and
        ps1 ... psn are error probabilities.
        Examples: {'rand': {'some_source': 0.1}} sets the random error rate
        for some_source to 0.1; {'rand': {''}}

    Attributes
    ----------
    prior_probs : dict[dict]
        A dictionary of prior systematic and random error probabilities for
        each knowledge source.
    subtype_probs: dict[dict]
        A dictionary of random error probabilities for knowledge sources.
        When a subtype random error probability is not specified, will just
        use the overall type prior in prior_probs. If None, will
        only use the priors for each rule.
    """
    def __init__(self, prior_probs=None, subtype_probs=None):
        self.prior_probs = default_probs
        if prior_probs:
            self.prior_probs.update(prior_probs)
        self.subtype_probs = subtype_probs

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
        if not use_reach_subtypes:
            logger.info('Belief engine could not import REACH subtypes, they '
                        'will be ignored.')
        self._check_prior_probs(statements)
        for st in statements:
            sources = [ev.source_api for ev in st.evidence]
            uniq_sources = numpy.unique(sources)
            syst_factors = {s: self.prior_probs['syst'][s]
                            for s in uniq_sources}
            rand_factors = {k: [] for k in uniq_sources}
            for ev in st.evidence:
                rand_factors[ev.source_api].append(
                        evidence_random_noise_prior(
                            ev,
                            self.prior_probs['rand'],
                            self.subtype_probs))

            neg_prob_prior = 1
            for s in uniq_sources:
                neg_prob_prior *= (syst_factors[s] +
                                   numpy.prod(rand_factors[s]))
            prob_prior = 1 - neg_prob_prior
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

    def _check_prior_probs(self, statements):
        """Check that we have probabilities  for sources in statements."""
        sources = set()
        for stmt in statements:
            sources |= set([ev.source_api for ev in stmt.evidence])
        for err_type in ('rand', 'syst'):
            for source in sources:
                if source not in self.prior_probs[err_type]:
                    msg = 'BeliefEngine missing probability parameter' + \
                        ' for source: %s' % source
                    raise Exception(msg)


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


def sample_statements(stmts, seed=None):
    """Return statements sampled according to belief.

    Statements are sampled independently according to their
    belief scores. For instance, a Staement with a belief
    score of 0.7 will end up in the returned Statement list
    with probability 0.7.

    Parameters
    ----------
    stmts : list[indra.statements.Statement]
        A list of INDRA Statements to sample.
    seed : Optional[int]
        A seed for the random number generator used for sampling.

    Returns
    -------
    new_stmts : list[indra.statements.Statement]
        A list of INDRA Statements that were chosen by random sampling
        according to their respective belief scores.
    """
    if seed:
        numpy.random.seed(seed)
    new_stmts = []
    r = numpy.random.rand(len(stmts))
    for i, stmt in enumerate(stmts):
        if r[i] < stmt.belief:
            new_stmts.append(stmt)
    return new_stmts


def evidence_random_noise_prior(evidence, type_probs, subtype_probs):
    """Determines the random-noise prior probability for this evidence.

    If the evidence corresponds to a subtype, and that subtype has a curated
    prior noise probability, use that.

    Otherwise, gives the random-noise prior for the overall rule type.
    """
    (stype, subtype) = tag_evidence_subtype(evidence)
    # Get the subtype, if available

    # Return the subtype random noise prior, if available
    if subtype_probs is not None:
        if stype in subtype_probs:
            if subtype in subtype_probs[stype]:
                return subtype_probs[stype][subtype]

    # Fallback to just returning the overall evidence type random noise prior
    return type_probs[stype]


def tag_evidence_subtype(evidence):
    """Returns the type and subtype of an evidence object as a string,
    typically the extraction rule or database from which the statement
    was generated.

    For biopax, this is just the database name.

    Parameters
    ----------
    statement: indra.statements.Evidence
        The statement which we wish to subtype

    Returns
    -------
    types: tuple
        A tuple with (type, subtype), both strings
        Returns (type, None) if the type of statement is not yet handled in
        this function.
    """

    source_api = evidence.source_api
    annotations = evidence.annotations

    if source_api == 'biopax':
        subtype = annotations.get('source_sub_id')
    elif source_api == 'reach':
        if 'found_by' in annotations:
            if use_reach_subtypes:
                subtype = determine_reach_subtype(annotations['found_by'])
            else:
                subtype = None
        else:
            logger.warning('Could not find found_by attribute in reach '
                           'statement annoations')
            subtype = None
    elif source_api == 'geneways':
        subtype = annotations['actiontype']
    else:
        subtype = None

    return (source_api, subtype)
