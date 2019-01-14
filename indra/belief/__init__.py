from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import json
import numpy
import logging
import networkx
from os import path, pardir
from collections import namedtuple


try:
    from indra.sources.reach.processor import determine_reach_subtype
    use_reach_subtypes = True
except ImportError:
    use_reach_subtypes = False

logger = logging.getLogger(__name__)


THIS_DIR = path.dirname(path.abspath(__file__))


def load_default_probs():
    json_path = path.join(THIS_DIR, pardir, 'resources',
                          'default_belief_probs.json') 
    with open(json_path, 'r') as f:
        prior_probs = json.load(f)
    return prior_probs


class BeliefScorer(object):
    """Base class for a belief engine scorer, which computes the
    prior probability of a statement given its type and evidence.

    To use with the belief engine, make a subclass with methods implemented.
    """
    def score_statement(self, st, extra_evidence=None):
        """Computes the prior belief probability for an INDRA Statement.

        The Statement is assumed to be de-duplicated. In other words,
        the Statement is assumed to have
        a list of Evidence objects that supports it. The prior probability of
        the Statement is calculated based on the number of Evidences it has
        and their sources.

        Parameters
        ----------
        st : indra.statements.Statement
            An INDRA Statements whose belief scores are to
            be calculated.
        extra_evidence : list[indra.statements.Evidence]
            A list of Evidences that are supporting the Statement (that aren't
            already included in the Statement's own evidence list.

        Returns
        -------
        belief_score : float
            The computed prior probability for the statement
        """
        raise NotImplementedError('Need to subclass BeliefScorer and '
                                  'implement methods.')

    def check_prior_probs(self, statements):
        """Make sure the scorer has all the information needed to compute
        belief scores of each statement in the provided list, and raises an
        exception otherwise.

        Parameters
        ----------
        statements : list<indra.statements.Statement>
            List of statements to check
        """
        raise NotImplementedError('Need to subclass BeliefScorer and '
                                  'implement methods.')


class SimpleScorer(BeliefScorer):
    """Computes the prior probability of a statement given its type and
    evidence.

    Parameters
    ----------
    prior_probs : dict[dict]
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
    subtype_probs : dict[dict]
        A dictionary of random error probabilities for knowledge sources.
        When a subtype random error probability is not specified, will just
        use the overall type prior in prior_probs. If None, will
        only use the priors for each rule.
    """
    def __init__(self, prior_probs=None, subtype_probs=None):
        self.prior_probs = load_default_probs()
        self.subtype_probs = {}
        self._update_probs(prior_probs, subtype_probs)

    def update_probs(prior_probs=None, subtype_probs=None):
        if prior_probs:
            for key in ('rand', 'syst'):
                self.prior_probs[key].update(prior_probs.get(key, {}))
        for err_type, source_dict in self.prior_probs.items():
            logger.debug("Prior probabilities for %s errors: %s"
                         % (err_type, source_dict))
        self.subtype_probs = subtype_probs

    def score_evidence_list(self, evidences):
        def _score(evidences):
            if not evidences:
                return 0
            # Collect all unique sources
            sources = [ev.source_api for ev in evidences]
            uniq_sources = numpy.unique(sources)
            # Calculate the systematic error factors given unique sources
            syst_factors = {s: self.prior_probs['syst'][s]
                            for s in uniq_sources}
            # Calculate the radom error factors for each source
            rand_factors = {k: [] for k in uniq_sources}
            for ev in evidences:
                rand_factors[ev.source_api].append(
                    evidence_random_noise_prior(
                        ev,
                        self.prior_probs['rand'],
                        self.subtype_probs))
            # The probability of incorrectness is the product of the
            # source-specific probabilities
            neg_prob_prior = 1
            for s in uniq_sources:
                neg_prob_prior *= (syst_factors[s] +
                                   numpy.prod(rand_factors[s]))
            # Finally, the probability of correctness is one minus incorrect
            prob_prior = 1 - neg_prob_prior
            return prob_prior
        pos_evidence = [ev for ev in evidences if
                        not ev.epistemics.get('negated')]
        neg_evidence = [ev for ev in evidences if
                        ev.epistemics.get('negated')]
        pp = _score(pos_evidence)
        np = _score(neg_evidence)
        # The basic assumption is that the positive and negative evidence
        # can't simultaneously be correct.
        # There are two cases to consider. (1) If the positive evidence is
        # incorrect then there is no Statement and the belief should be 0,
        # irrespective of the negative evidence.
        # (2) If the positive evidence is correct and the negative evidence
        # is incorrect.
        # This amounts to the following formula:
        # 0 * (1-pp) + 1 * (pp * (1-np)) which we simplify below
        score = pp * (1 - np)
        return score

    def score_statement(self, st, extra_evidence=None):
        """Computes the prior belief probability for an INDRA Statement.

        The Statement is assumed to be de-duplicated. In other words,
        the Statement is assumed to have
        a list of Evidence objects that supports it. The prior probability of
        the Statement is calculated based on the number of Evidences it has
        and their sources.

        Parameters
        ----------
        st : indra.statements.Statement
            An INDRA Statements whose belief scores are to
            be calculated.
        extra_evidence : list[indra.statements.Evidence]
            A list of Evidences that are supporting the Statement (that aren't
            already included in the Statement's own evidence list.

        Returns
        -------
        belief_score : float
            The computed prior probability for the statement
        """
        if extra_evidence is None:
            extra_evidence = []
        all_evidence = st.evidence + extra_evidence
        return self.score_evidence_list(all_evidence)

    def check_prior_probs(self, statements):
        """Throw Exception if BeliefEngine parameter is missing.

        Make sure the scorer has all the information needed to compute
        belief scores of each statement in the provided list, and raises an
        exception otherwise.

        Parameters
        ----------
        statements : list[indra.statements.Statement]
            List of statements to check
        """
        sources = set()
        for stmt in statements:
            sources |= set([ev.source_api for ev in stmt.evidence])
        for err_type in ('rand', 'syst'):
            for source in sources:
                if source not in self.prior_probs[err_type]:
                    msg = 'BeliefEngine missing probability parameter' + \
                        ' for source: %s' % source
                    raise Exception(msg)


default_scorer = SimpleScorer()


class BeliefEngine(object):
    """Assigns beliefs to INDRA Statements based on supporting evidence.

    Attributes
    ----------
    scorer : BeliefScorer
        A BeliefScorer object that computes the prior probability of a
        statement given its its statment type and evidence.
        Must implement the `score_statement` method which takes
        Statements and computes the belief score of a statement, and the
        `check_prior_probs` method which takes a list of INDRA Statements and
        verifies that the scorer has all the information it needs to score
        every statement in the list, and raises an exception if not.
    """
    def __init__(self, scorer=None):
        if scorer is None:
            scorer = default_scorer
        assert(isinstance(scorer, BeliefScorer))
        self.scorer = scorer

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
        self.scorer.check_prior_probs(statements)
        for st in statements:
            st.belief = self.scorer.score_statement(st)

    def set_hierarchy_probs(self, statements):
        """Sets hierarchical belief probabilities for INDRA Statements.

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
            """Return a DiGraph based on matches keys and Statement supports"""
            g = networkx.DiGraph()
            for st1 in stmts:
                g.add_node(st1.matches_key(), stmt=st1)
                for st2 in st1.supported_by:
                    g.add_node(st2.matches_key(), stmt=st2)
                    g.add_edge(st2.matches_key(), st1.matches_key())
            return g

        def get_ranked_stmts(g):
            """Return a topological sort of statement matches keys from a graph.
            """
            node_ranks = networkx.algorithms.dag.topological_sort(g)
            node_ranks = reversed(list(node_ranks))
            stmts = [g.node[n]['stmt'] for n in node_ranks]
            return stmts

        def assert_no_cycle(g):
            """If the graph has cycles, throws AssertionError."""
            try:
                cyc = networkx.algorithms.cycles.find_cycle(g)
            except networkx.exception.NetworkXNoCycle:
                return
            msg = 'Cycle found in hierarchy graph: %s' % cyc
            assert False, msg

        g = build_hierarchy_graph(statements)
        assert_no_cycle(g)
        ranked_stmts = get_ranked_stmts(g)
        for st in ranked_stmts:
            bps = _get_belief_package(st)
            supporting_evidences = []
            # NOTE: the last belief package in the list is this statement's own
            for bp in bps[:-1]:
                # Iterate over all the parent evidences and add only
                # non-negated ones
                for ev in bp.evidences:
                    if not ev.epistemics.get('negated'):
                        supporting_evidences.append(ev)
            # Now add the Statement's own evidence
            # Now score all the evidences
            belief = self.scorer.score_statement(st, supporting_evidences)
            st.belief = belief

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


BeliefPackage = namedtuple('BeliefPackage', 'statement_key evidences')


def _get_belief_package(stmt):
    """Return the belief packages of a given statement recursively."""
    # This list will contain the belief packages for the given statement
    belief_packages = []
    # Iterate over all the support parents
    for st in stmt.supports:
        # Recursively get all the belief packages of the parent
        parent_packages = _get_belief_package(st)
        package_stmt_keys = [pkg.statement_key for pkg in belief_packages]
        for package in parent_packages:
            # Only add this belief package if it hasn't already been added
            if package.statement_key not in package_stmt_keys:
                belief_packages.append(package)
    # Now make the Statement's own belief package and append it to the list
    belief_package = BeliefPackage(stmt.matches_key(), stmt.evidence)
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
    elif source_api in ('reach', 'eidos'):
        if 'found_by' in annotations:
            if source_api == 'reach' and use_reach_subtypes:
                subtype = determine_reach_subtype(annotations['found_by'])
            elif source_api == 'eidos':
                subtype = annotations['found_by']
            else:
                subtype = None
        else:
            logger.debug('Could not find found_by attribute in reach '
                         'statement annoations')
            subtype = None
    elif source_api == 'geneways':
        subtype = annotations['actiontype']
    else:
        subtype = None

    return (source_api, subtype)
