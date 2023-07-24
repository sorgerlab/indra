import copy
import json
import tqdm
import numpy
import logging
import networkx
from os import path, pardir
from typing import List, Optional, Dict, Callable, Tuple, Sequence, Iterable
from indra.mechlinker import LinkedStatement
from indra.statements import Evidence, Statement


logger = logging.getLogger(__name__)


THIS_DIR = path.dirname(path.abspath(__file__))


def load_default_probs() -> Dict[str, Dict[str, float]]:
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
    def score_statements(
        self,
        statements: Sequence[Statement],
        extra_evidence: Optional[List[List[Evidence]]] = None,
    ) -> Sequence[float]:
        """Computes belief probabilities for a list of INDRA Statements.

        The Statements are assumed to be de-duplicated. In other words, each
        Statement is assumed to have a list of Evidence objects that supports
        it. The probability of correctness of the Statement is generally
        calculated based on the number of Evidences it has, their sources, and
        other features depending on the subclass implementation.

        Parameters
        ----------
        statements :
            INDRA Statements whose belief scores are to be calculated.
        extra_evidence :
            A list corresponding to the given list of statements, where
            each entry is a list of Evidence objects providing additional
            support for the corresponding statement (i.e., Evidences that
            aren't already included in the Statement's own evidence list).

        Returns
        -------
        :
            The computed prior probabilities for each statement.
        """
        raise NotImplementedError('Need to subclass BeliefScorer and '
                                  'implement methods.')

    def score_statement(
        self,
        statement: Statement,
        extra_evidence: Optional[List[Evidence]] = None,
    ) -> float:
        """Score a single statement by passing arguments to `score_statements`.
        """
        extra_evidence_wrap = None if extra_evidence is None \
                                   else [extra_evidence]
        return self.score_statements([statement], extra_evidence_wrap)[0]

    def check_prior_probs(
        self,
        statements: Sequence[Statement],
    ) -> None:
        """Make sure the scorer has all the information needed to compute
        belief scores of each statement in the provided list, and raises an
        exception otherwise.

        Parameters
        ----------
        statements :
            List of statements to check
        """
        raise NotImplementedError('Need to subclass BeliefScorer and '
                                  'implement methods.')


class SimpleScorer(BeliefScorer):
    """Computes the prior probability of a statement given its type and
    evidence.

    Parameters
    ----------
    prior_probs :
        A dictionary of prior probabilities used to override/extend the default
        ones. There are two types of prior probabilities: rand and syst,
        corresponding to random error and systematic error rate for each
        knowledge source. The prior_probs dictionary has the general structure
        {'rand': {'s1': pr1, ..., 'sn': prn}, 'syst': {'s1': ps1, ..., 'sn':
        psn}} where 's1' ... 'sn' are names of input sources and pr1 ... prn
        and ps1 ... psn are error probabilities.  Examples: {'rand':
        {'some_source': 0.1}} sets the random error rate for some_source to
        0.1; {'rand': {''}}
    subtype_probs :
        A dictionary of random error probabilities for knowledge sources.
        When a subtype random error probability is not specified, will just
        use the overall type prior in prior_probs. If None, will
        only use the priors for each rule.
    """
    def __init__(
        self,
        prior_probs: Optional[Dict[str, Dict[str, float]]] = None,
        subtype_probs: Optional[Dict[str, Dict[str, float]]] = None,
    ):
        self.prior_probs = load_default_probs()
        self.subtype_probs: Optional[Dict[str, Dict[str, float]]] = {}
        self.update_probs(prior_probs, subtype_probs)

    def update_probs(
        self,
        prior_probs: Optional[Dict[str, Dict[str, float]]] = None,
        subtype_probs: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> None:
        """Update Scorer's prior probabilities with the given dictionaries."""
        if prior_probs:
            for key in ('rand', 'syst'):
                self.prior_probs[key].update(prior_probs.get(key, {}))
        for err_type, source_dict in self.prior_probs.items():
            logger.debug("Prior probabilities for %s errors: %s"
                         % (err_type, source_dict))
        self.subtype_probs = subtype_probs

    def score_evidence_list(
        self,
        evidences: List[Evidence],
    ) -> float:
        """Return belief score given a list of supporting evidences.

        Parameters
        ----------
        evidences :
            List of evidences to use for calculating a statement's belief.

        Returns
        -------
        :
            Belief value based on the evidences.
        """
        def _score(evidences):
            if not evidences:
                return 0
            # Collect all unique sources
            sources = [ev.source_api for ev in evidences]
            uniq_sources = numpy.unique(sources)
            # Calculate the systematic error factors given unique sources
            syst_factors = {s: self.prior_probs['syst'][s]
                            for s in uniq_sources}
            # Calculate the random error factors for each source
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
        # Split evidence into positive and negative and score
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

    def score_statements(
        self,
        statements: Sequence[Statement],
        extra_evidence: Optional[List[List[Evidence]]] = None,
    ) -> List[float]:
        """Computes belief probabilities for a list of INDRA Statements.

        The Statements are assumed to be de-duplicated. In other words, each
        Statement is assumed to have a list of Evidence objects that supports
        it. The probability of correctness of the Statement is generally
        calculated based on the number of Evidences it has, their sources, and
        other features depending on the subclass implementation.

        Parameters
        ----------
        statements :
            INDRA Statements whose belief scores are to be calculated.
        extra_evidence :
            A list corresponding to the given list of statements, where
            each entry is a list of Evidence objects providing additional
            support for the corresponding statement (i.e., Evidences that
            aren't already included in the Statement's own evidence list).

        Returns
        -------
        :
            The computed prior probabilities for each statement.
        """
        # Check our list of extra evidences
        check_extra_evidence(extra_evidence, len(statements))
        # Get beliefs for each statement
        beliefs = []
        for ix, stmt in enumerate(statements):
            all_evidence = get_stmt_evidence(stmt, ix, extra_evidence)
            beliefs.append(self.score_evidence_list(all_evidence))
        return beliefs

    def check_prior_probs(
        self,
        statements: Sequence[Statement],
    ) -> None:
        """Throw Exception if BeliefEngine parameter is missing.

        Make sure the scorer has all the information needed to compute
        belief scores of each statement in the provided list, and raises an
        exception otherwise.

        Parameters
        ----------
        statements :
            List of statements to check.
        """
        sources = set()
        for stmt in statements:
            sources |= set([ev.source_api for ev in stmt.evidence])
        return self._check_sources(sources)

    def _check_sources(
        self,
        sources: Iterable[str],
    ) -> None:
        """Make sure all sources have entries for prior parameters."""
        for err_type in ('rand', 'syst'):
            for source in sources:
                if source not in self.prior_probs[err_type]:
                    msg = 'BeliefEngine missing probability parameter' + \
                        ' for source: %s' % source
                    raise Exception(msg)


class BayesianScorer(SimpleScorer):
    """This is a belief scorer which assumes a Beta prior and a set of prior
    counts of correct and incorrect instances for a given source. It exposes
    and interface to take additional counts and update its probability
    parameters which can then be used to calculate beliefs on a set of
    Statements.

    Parameters
    ----------
    prior_counts :
        A dictionary of counts of the form [pos, neg] for
        each source.
    subtype_counts :
        A dict of dicts of counts of the form [pos, neg] for each subtype
        within a source.
    """
    def __init__(
        self,
        prior_counts: Dict[str, List[int]],
        subtype_counts: Dict[str, Dict[str, List[int]]],
    ):
        self.prior_probs = load_default_probs()
        self.subtype_probs = {}
        self.prior_counts = copy.deepcopy(prior_counts) if prior_counts else {}
        self.subtype_counts = copy.deepcopy(subtype_counts) if subtype_counts \
            else {}
        # Set the probability estimates based on the counts
        self.update_probs()

    def update_probs(self):
        """Update the internal probability values given the counts."""
        # We deal with the prior probsfirst
        # This is a fixed assumed value for systematic error
        syst_error = 0.05
        prior_probs = {'syst': {}, 'rand': {}}
        for source, (p, n) in self.prior_counts.items():
            # Skip if there are no actual counts
            if n + p == 0:
                continue
            prior_probs['syst'][source] = syst_error
            prior_probs['rand'][source] = \
                1 - min((float(p) / (n + p), 1-syst_error)) - syst_error
        # Next we deal with subtype probs based on counts
        subtype_probs = {}
        for source, entry in self.subtype_counts.items():
            for rule, (p, n) in entry.items():
                # Skip if there are no actual counts
                if n + p == 0:
                    continue
                if source not in subtype_probs:
                    subtype_probs[source] = {}
                subtype_probs[source][rule] = \
                    1 - min((float(p) / (n + p), 1-syst_error)) - syst_error
        # Finally we propagate this into the full probability
        # data structures of the parent class
        super(BayesianScorer, self).update_probs(prior_probs, subtype_probs)

    def update_counts(
        self,
        prior_counts: Dict[str, List[int]],
        subtype_counts: Dict[str, Dict[str, List[int]]],
    ) -> None:
        """Update the internal counts based on given new counts.

        Parameters
        ----------
        prior_counts :
            A dictionary of counts of the form [pos, neg] for
            each source.
        subtype_counts :
            A dict of dicts of counts of the form [pos, neg] for each subtype
            within a source.
        """
        for source, (pos, neg) in prior_counts.items():
            if source not in self.prior_counts:
                self.prior_counts[source] = [0, 0]
            self.prior_counts[source][0] += pos
            self.prior_counts[source][1] += neg
        for source, subtype_dict in subtype_counts.items():
            if source not in self.subtype_counts:
                self.subtype_counts[source] = {}
            for subtype, (pos, neg) in subtype_dict.items():
                if subtype not in self.subtype_counts[source]:
                    self.subtype_counts[source][subtype] = [0, 0]
                self.subtype_counts[source][subtype][0] += pos
                self.subtype_counts[source][subtype][1] += neg
        self.update_probs()


default_scorer = SimpleScorer()


class BeliefEngine(object):
    """Assigns beliefs to INDRA Statements based on supporting evidence.

    Parameters
    ----------
    scorer :
        A BeliefScorer object that computes the prior probability of a
        statement given its its statment type, evidence, or other features.
        Must implement the `score_statements` method which takes Statements and
        computes the belief score of a statement, and the `check_prior_probs`
        method which takes a list of INDRA Statements and verifies that the
        scorer has all the information it needs to score every statement in the
        list, and raises an exception if not.
    matches_fun :
        A function handle for a custom matches key if a non-default one is
        used. Default is None.
    refinements_graph :
        A graph whose nodes are statement hashes, and edges point from a more
        specific to a less specific statement representing a refinement. If not
        given, a new graph is constructed here.
    """
    def __init__(
        self,
        scorer: Optional[BeliefScorer] = None,
        matches_fun: Optional[Callable[[Statement], str]] = None,
        refinements_graph: Optional[networkx.DiGraph] = None,
    ):
        if scorer is None:
            scorer = default_scorer
        assert isinstance(scorer, BeliefScorer)
        self.scorer = scorer

        self.matches_fun = matches_fun if matches_fun else \
            lambda stmt: stmt.matches_key()

        self.refinements_graph = refinements_graph

    def set_prior_probs(
        self,
        statements: Sequence[Statement],
    ) -> None:
        """Sets the prior belief probabilities for a list of INDRA Statements.

        The Statements are assumed to be de-duplicated. In other words,
        each Statement in the list passed to this function is assumed to have
        a list of Evidence objects that support it. The prior probability of
        each Statement is calculated based on the number of Evidences it has
        and their sources.

        Parameters
        ----------
        statements :
            A list of INDRA Statements whose belief scores are to
            be calculated. Each Statement object's belief attribute is updated
            by this function.
        """
        self.scorer.check_prior_probs(statements)
        prior_probs = self.scorer.score_statements(statements)
        assert len(prior_probs) == len(statements), \
                "prior_probs and statements should be the same length"
        for ix, st in enumerate(statements):
            st.belief = prior_probs[ix]

    def set_hierarchy_probs(
        self,
        statements: Sequence[Statement],
    ) -> None:
        """Sets hierarchical belief probabilities for INDRA Statements.

        The Statements are assumed to be in a hierarchical relation graph with
        the supports and supported_by attribute of each Statement object having
        been set. The hierarchical belief probability of each Statement is
        calculated based the accumulated evidence from both itself and its more
        specific statements in the hierarchy graph.

        Parameters
        ----------
        statements :
            A list of INDRA Statements whose belief scores are to
            be calculated. Each Statement object's belief attribute is updated
            by this function.
        """
        beliefs = self.get_hierarchy_probs(statements)
        for stmt in statements:
            sh = stmt.get_hash(matches_fun=self.matches_fun)
            stmt.belief = beliefs[sh]

    def get_hierarchy_probs(
        self,
        statements: Sequence[Statement],
    ) -> Dict[int, float]:
        """Gets hierarchical belief probabilities for INDRA Statements.

        Parameters
        ----------
        statements :
            A list of INDRA Statements whose belief scores are to
            be calculated. Each Statement object's belief attribute is updated
            by this function.

        Returns
        -------
        :
            A dictionary mapping statement hashes to corresponding belief
            scores. Hashes are calculated using the instance's `self.matches_fun`.
        """
        # We only re-build the refinements graph if one wasn't provided
        # as an argument
        if self.refinements_graph is None:
            # Build the graph for the given set of statements
            self.refinements_graph = build_refinements_graph(statements,
                                                   matches_fun=self.matches_fun)
        # Get the evidences from the more specific (supports) statements
        all_extra_evs = get_ev_for_stmts_from_supports(statements,
                                                   self.refinements_graph)
        return self._hierarchy_probs_from_evidences(statements, all_extra_evs)

    def get_hierarchy_probs_from_hashes(
        self,
        statements: Sequence[Statement],
        refiners_list: List[List[int]],
    ) -> Dict[int, float]:
        """Return the full belief of a statement with refiners given as hashes.

        Parameters
        ----------
        statements :
            Statements to calculate beliefs for.
        refiners_list :
            A list corresponding to the list of statements, where each entry
            is a list of statement hashes for the statements that are
            refinements (i.e., more specific versions) of the corresponding
            statement in the statements list. If there are no refiner
            statements the entry should be an empty list.

        Returns
        -------
        :
            A dictionary mapping statement hashes to corresponding belief
            scores.
        """
        if self.refinements_graph is None:
            raise ValueError("refinements_graph not initialized.")
        # Get the evidences from the more specific (supports) statements
        all_extra_evs = get_ev_for_stmts_from_hashes(statements,
                                                     refiners_list,
                                                     self.refinements_graph)
        # Return beliefs using all the evidences
        return self._hierarchy_probs_from_evidences(statements, all_extra_evs)

    def _hierarchy_probs_from_evidences(
        self,
        statements: Sequence[Statement],
        extra_evidence: List[List[Evidence]],
    ) -> Dict[int, float]:
        """Use the Scorer to get stmt beliefs with supports evidences."""
        # Get the list of beliefs matching the statements we passed in
        beliefs = self.scorer.score_statements(statements, extra_evidence)
        # Convert to a dict of beliefs keyed by hash and return
        hashes = [s.get_hash(matches_fun=self.matches_fun) for s in statements]
        beliefs_by_hash = dict(zip(hashes, beliefs))
        return beliefs_by_hash

    def set_linked_probs(
        self,
        linked_statements: List[LinkedStatement],
    ) -> None:
        """Sets the belief probabilities for a list of linked INDRA Statements.

        The list of LinkedStatement objects is assumed to come from the
        MechanismLinker. The belief probability of the inferred Statement is
        assigned the joint probability of its source Statements.

        Parameters
        ----------
        linked_statements :
            A list of INDRA LinkedStatements whose belief scores are to
            be calculated. The belief attribute of the inferred Statement in
            the LinkedStatement object is updated by this function.
        """
        for st in linked_statements:
            source_probs = [s.belief for s in st.source_stmts]
            st.inferred_stmt.belief = numpy.prod(source_probs)


def get_ev_for_stmts_from_supports(
    statements: Sequence[Statement],
    refinements_graph: Optional[networkx.DiGraph] = None,
    matches_fun: Optional[Callable[[Statement], str]] = None,
) -> List[List[Evidence]]:
    """
    Collect evidence from the more specific statements of a list of statements.

    Parameters
    ----------
    statements :
        A list of Statements with `supports` Statements.
    refinements_graph :
        A networkx graph whose nodes are statement hashes carrying a stmt
        attribute with the actual statement object. Edges point from less
        detailed to more detailed statements (i.e., from a statement to another
        statement that refines it). If not provided, the graph is generated
        from `statements` using the function
        :py:func:`build_refinements_graph`.
    matches_fun :
        An optional function to calculate the matches key and hash of a
        given statement. If not provided, the default matches function is
        used. Default: None.

    Returns
    -------
    :
        A list corresponding to the given list of statements, where each entry
        is a list of Evidence objects providing additional support for the
        corresponding statement (i.e., Evidences that aren't already included
        in the Statement's own evidence list).
    """
    # If the refinements_graph was not given, build it for this set of
    # statements
    if refinements_graph is None:
        # Build the graph for the given set of statements
        refinements_graph = build_refinements_graph(statements,
                                                    matches_fun=matches_fun)
    try:
        assert_no_cycle(refinements_graph)
    except AssertionError as err:
        cycle_file = path.join(THIS_DIR, 'refinement_cycles')
        logger.debug(f'Cycles found. Saving them to {cycle_file}, '
                     f'then raising error')
        find_cycles(g=refinements_graph, fpath=cycle_file)
        raise err

    # Use graph to collect corresponding lists of refiners for the statements
    refiners_list = []
    for stmt in statements:
        stmt_hash = stmt.get_hash(matches_fun=matches_fun)
        # Get the refiners/more specific stmts, if any (the edges in the
        # graph point from general to specific):
        refiners = list(networkx.descendants(refinements_graph, stmt_hash))
        refiners_list.append(refiners)
    # Use the refiners hashes to get the evidences
    return get_ev_for_stmts_from_hashes(statements, refiners_list,
                                        refinements_graph)


def get_ev_for_stmts_from_hashes(
    statements: Sequence[Statement],
    refiners_list: List[List[int]],
    refinements_graph: networkx.DiGraph,
) -> List[List[Evidence]]:
    """
    Collect evidence from the more specific statements of a list of statements.

    Similar to :py:func:`get_ev_for_stmts_from_supports`, but the more specific
    statements are specified explicitly by the hashes in `refiners_list` rather
    than obtained from the `supports` attribute of each statement. In addition,
    the `refinements_graph` argument is expected to have been pre-calculated
    (using the same matches key function used to generate the hashes in
    `refiners_list`) and hence is not optional.

    Parameters
    ----------
    statements :
        A list of Statements with `supports` Statements.
    refiners_list :
        A list corresponding to the list of statements, where each entry
        is a list of statement hashes for the statements that are
        refinements (i.e., more specific versions) of the corresponding
        statement in the statements list. If there are no refiner
        statements for a statement the entry should be an empty list.
    refinements_graph :
        A networkx graph whose nodes are statement hashes carrying a stmt
        attribute with the actual statement object. Edges point from less
        detailed to more detailed statements (i.e., from a statement to another
        statement that refines it).

    Returns
    -------
    :
        A list corresponding to the given list of statements, where each entry
        is a list of Evidence objects providing additional support for the
        corresponding statement (i.e., Evidences that aren't already included
        in the Statement's own evidence list).
    """
    all_extra_evs = []
    for stmt, refiners in zip(statements, refiners_list):
        # Collect evidence from all refiners while excluding any negated
        # evidence. Negated evidence for a more specific statement
        # isn't currently considered as counting against the believability
        # of a more general statement.
        extra_ev_for_stmt = list(set(
            ev
            for supp in refiners
            for ev in refinements_graph.nodes[supp]['stmt'].evidence
            if not ev.epistemics.get('negated')
        ))
        # Add the extra evidences for the statement to the full list
        all_extra_evs.append(extra_ev_for_stmt)
    return all_extra_evs


def sample_statements(
    stmts: Sequence[Statement],
    seed: Optional[int] = None,
) -> List[Statement]:
    """Return statements sampled according to belief.

    Statements are sampled independently according to their belief scores. For
    instance, a Statement with a belief score of 0.7 will end up in the returned
    Statement list with probability 0.7.

    Parameters
    ----------
    stmts :
        A list of INDRA Statements to sample.
    seed :
        A seed for the random number generator used for sampling.

    Returns
    -------
    :
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


def evidence_random_noise_prior(
    evidence: Evidence,
    type_probs: Dict[str, float],
    subtype_probs: Optional[Dict[str, Dict[str, float]]],
) -> float:
    """Gets the random-noise prior probability for this evidence.

    If the evidence corresponds to a subtype, and that subtype has a curated
    prior noise probability, use that.

    Otherwise, gives the random-noise prior for the overall rule type.
    """
    # Get the subtype, if available
    (stype, subtype) = tag_evidence_subtype(evidence)

    # Return the subtype random noise prior, if available
    if subtype_probs is not None:
        if stype in subtype_probs:
            if subtype in subtype_probs[stype]:
                return subtype_probs[stype][subtype]

    # Fallback to just returning the overall evidence type random noise prior
    return type_probs[stype]


def tag_evidence_subtype(
    evidence: Evidence,
) -> Tuple[str, Optional[str]]:
    """Returns the type and subtype of an evidence object as a string,
    typically the extraction rule or database from which the statement
    was generated.

    For biopax, this is just the database name.

    Parameters
    ----------
    statement:
        The statement which we wish to subtype

    Returns
    -------
    :
        A tuple with (type, subtype), both strings. Returns (type, None) if the
        type of statement is not yet handled in this function.
    """
    source_api = evidence.source_api
    annotations = evidence.annotations

    if source_api == 'biopax':
        subtype = annotations.get('source_sub_id')
    elif source_api in ('reach', 'eidos'):
        if 'found_by' in annotations:
            from indra.sources.reach.processor import determine_reach_subtype
            if source_api == 'reach':
                subtype = determine_reach_subtype(annotations['found_by'])
            elif source_api == 'eidos':
                subtype = annotations['found_by']
            else:
                subtype = None
        else:
            logger.debug('Could not find found_by attribute in reach '
                         'statement annotations')
            subtype = None
    elif source_api == 'geneways':
        subtype = annotations['actiontype']
    else:
        subtype = None

    return (source_api, subtype)


def build_refinements_graph(
    statements: Sequence[Statement],
    matches_fun: Optional[Callable[[Statement], str]] = None,
) -> networkx.DiGraph:
    """Return a DiGraph based on matches hashes and Statement refinements.

    Parameters
    ----------
    statements :
        A list of Statements with `supports` Statements, used to generate the
        refinements graph.
    matches_fun :
        An optional function to calculate the matches key and hash of a
        given statement. Default: None

    Returns
    -------
    :
        A networkx graph whose nodes are statement hashes carrying a stmt
        attribute with the actual statement object. Edges point from less
        detailed to more detailed statements (i.e., from a statement to another
        statement that refines it).
    """
    logger.debug('Building refinements graph')
    g = networkx.DiGraph()
    for st1 in statements:
        sh1 = st1.get_hash(matches_fun=matches_fun)
        g.add_node(sh1, stmt=st1)
        for st2 in st1.supports:
            sh2 = st2.get_hash(matches_fun=matches_fun)
            g.add_node(sh2, stmt=st2)
            g.add_edge(sh1, sh2)
    logger.debug('Finished building refinements graph')
    return g


def extend_refinements_graph(
    g: networkx.DiGraph,
    stmt: Statement,
    less_specifics: List[int],
    matches_fun: Optional[Callable[[Statement], str]] = None,
) -> networkx.DiGraph:
    """Extend refinements graph with a new statement and its refinements.

    Parameters
    ----------
    g :
        A refinements graph to be extended.
    stmt :
        The statement to be added to the refinements graph.
    less_specifics :
        A list of statement hashes of statements that are refined
        by this statement (i.e., are less specific versions of it).
    matches_fun :
        An optional function to calculate the matches key and hash of a
        given statement. Default: None
    """
    sh = stmt.get_hash(matches_fun=matches_fun)
    g.add_node(sh, stmt=stmt)
    for less_spec in less_specifics:
        g.add_edge(less_spec, sh)
    return g


def assert_no_cycle(
    g: networkx.DiGraph
) -> None:
    """If the graph has cycles, throws AssertionError.

    This can be used to make sure that a refinements graph is a DAG.

    Parameters
    ----------
    g :
        A refinements graph.
    """
    logger.debug('Looking for cycles in belief graph')
    try:
        cyc = networkx.algorithms.cycles.find_cycle(g)
    except networkx.exception.NetworkXNoCycle:
        return
    msg = 'Cycle found in hierarchy graph: %s' % cyc
    assert False, msg


def find_cycles(
    g: networkx.DiGraph,
    fpath: str,
    upload_to_s3: bool = True
) -> None:
    from datetime import datetime
    logger.debug('Looking for cycles')

    # Create cycle generator, see:
    # https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.cycles.simple_cycles.html#networkx.algorithms.cycles.simple_cycles
    cyc_gen = networkx.algorithms.cycles.simple_cycles(g)
    cycles = []
    for cyc in tqdm.tqdm(cyc_gen):
        cycles.append(cyc)

    if cycles:
        dt = datetime.now()
        with open(fpath, 'w') as fp:
            for c in cycles:
                # Each cycle is a list of lists
                cs = ','.join([str(sc) for sc in c])
                fp.write(f'{cs}\n')
        logger.debug(f'Cycles written to {fpath}')
        if upload_to_s3:
            from indra.util.aws import get_s3_client
            s3 = get_s3_client(unsigned=False)
            fname = f'indra-db/dumps/refinement_cycles_' \
                    f'{dt.strftime("%Y%m%d%H%M%S")}'
            s3.upload_file(fpath, 'bigmech', fname)
            logger.debug(f'Cycles uploaded to {fname}')
    else:
        logger.debug('No cycles were found')


def get_ranked_stmts(g):
    """Return a topological sort of statements from a graph."""
    logger.debug('Getting ranked statements')
    node_ranks = networkx.algorithms.dag.topological_sort(g)
    node_ranks = reversed(list(node_ranks))
    stmts = [g.nodes[n]['stmt'] for n in node_ranks]
    return stmts


def check_extra_evidence(
    extra_evidence: Optional[List[List[Evidence]]],
    num_stmts: int,
) -> None:
    """Check whether extra evidence list has correct length/contents.

    Raises ValueError if the extra_evidence list does not match the length
    num_stmts, or if it contains items other than empty lists or lists
    of Evidence objects.

    Parameters
    ----------
    extra_evidence :
        A list of length num_stmts where each entry is a list of Evidence
        objects, or None. If extra_evidence is None, the function returns
        without raising an error.
    num_stmts :
        An integer giving the required length of the extra_evidence list
        (which should correspond to a list of statements)
    """
    # If given, check the extra_evidence list
    if extra_evidence is not None:
        if num_stmts != len(extra_evidence):
            raise ValueError("extra_evidence must be a list of the same "
                             "length as stmts.")
        for entry in extra_evidence:
            if not (isinstance(entry, list)):
                raise ValueError("extra_evidence must be a list of lists.")


def get_stmt_evidence(
    stmt: Statement,
    ix: int,
    extra_evidence: Optional[List[List[Evidence]]],
) -> List[Evidence]:
    """Combine a statements' own evidence with any extra evidence provided."""
    stmt_ev = set(stmt.evidence)
    if extra_evidence is not None:
        extra_ev_for_stmt = extra_evidence[ix]
        stmt_ev.update(extra_ev_for_stmt)
    return list(stmt_ev)



