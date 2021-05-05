import copy
import json
import numpy
import logging
import networkx
from os import path, pardir


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
        self.update_probs(prior_probs, subtype_probs)

    def update_probs(self, prior_probs=None, subtype_probs=None):
        if prior_probs:
            for key in ('rand', 'syst'):
                self.prior_probs[key].update(prior_probs.get(key, {}))
        for err_type, source_dict in self.prior_probs.items():
            logger.debug("Prior probabilities for %s errors: %s"
                         % (err_type, source_dict))
        self.subtype_probs = subtype_probs

    def score_evidence_list(self, evidences):
        """Return belief score given a list of supporting evidences."""
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
        # We remove instance duplicates here
        all_evidence = set(st.evidence) | set(extra_evidence)
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


class BayesianScorer(SimpleScorer):
    """This is a belief scorer which assumes a Beta prior and a set of prior
    counts of correct and incorrect instances for a given source. It exposes
    and interface to take additional counts and update its probability
    parameters which can then be used to calculate beliefs on a set of
    Statements.

    Parameters
    ----------
    prior_counts : dict
        A dictionary of counts of the form [pos, neg] for
        each source.
    subtype_counts : dict
        A dictionary of counts of the form [pos, neg] for
        each subtype within a source.
    """
    def __init__(self, prior_counts, subtype_counts):
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

    def update_counts(self, prior_counts, subtype_counts):
        """Update the internal counts based on given new counts.

        Parameters
        ----------
        prior_counts : dict
            A dictionary of counts of the form [pos, neg] for
            each source.
        subtype_counts : dict
            A dictionary of counts of the form [pos, neg] for
            each subtype within a source.
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

    Attributes
    ----------
    scorer : Optional[BeliefScorer]
        A BeliefScorer object that computes the prior probability of a
        statement given its its statment type and evidence.
        Must implement the `score_statement` method which takes
        Statements and computes the belief score of a statement, and the
        `check_prior_probs` method which takes a list of INDRA Statements and
        verifies that the scorer has all the information it needs to score
        every statement in the list, and raises an exception if not.
    matches_fun : Optional[function]
        A function handle for a custom matches key if a non-deafult one is
        used. Default: None
    refinements_graph : Optional[networkx.DiGraph]
        A graph whose nodes are statement hashes, and edges point from
        a more specific to a less specific statement representing
        a refinement. If not given, a new graph is constructed here.
    """
    def __init__(self, scorer=None, matches_fun=None, refinements_graph=None):
        if scorer is None:
            scorer = default_scorer
        assert isinstance(scorer, BeliefScorer)
        self.scorer = scorer

        self.matches_fun = matches_fun if matches_fun else \
            lambda stmt: stmt.matches_key()

        self.refinements_graph = refinements_graph

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
        self.scorer.check_prior_probs(statements)
        prior_probs = self.scorer.score_statements(statements)
        assert len(prior_probs) == len(statements), \
                "prior_probs and statements should be the same length"
        for ix, st in enumerate(statements):
            st.belief = prior_probs[ix]


    def get_refinement_probs(self, statements, refiners_list=None):
        """Return the full belief of a statement given its refiners.

        Parameters
        ----------
        statements : list of indra.statements.Statement
            Statements to calculate beliefs for.
        refiners_list : list[list[int]]
            A list corresponding to the list of statements, where each entry
            is a list of statement hashes for the statements that are
            refinements (i.e., more specific versions) of the corresponding
            statement in the statements list. If there are no refiner
            statements the entry should be an empty list.

        Returns
        -------
        ****FIXME FIXME FIXME
        iterable of floats
            Belief scores for the list of statements.
        """
        all_extra_evs = []
        for stmt, refiners in zip(statements, refiners_list):
            extra_ev_for_stmt = set()
            for supp in refiners:
                extra_ev_for_stmt |= \
                    set(self.refinements_graph.nodes[supp]['stmt'].evidence)
            # Exclude any negated evidences
            extra_ev_for_stmt = list({ev for ev in extra_ev_for_stmt
                                         if not ev.epistemics.get('negated')})
            all_extra_evs.append(extra_ev_for_stmt)

        # TODO
        #beliefs = self.scorer.score_statements(statements, all_extra_evs)
        beliefs = self.scorer.score_statements(statements)
        hashes = [s.get_hash(self.matches_fun) for s in statements]
        beliefs_by_hash = dict(zip(hashes, beliefs))
        return beliefs_by_hash

    def set_hierarchy_probs(self, statements):
        """Sets hierarchical belief probabilities for INDRA Statements.

        The Statements are assumed to be in a hierarchical relation graph with
        the supports and supported_by attribute of each Statement object having
        been set.
        The hierarchical belief probability of each Statement is calculated
        based on its prior probability and the probabilities propagated from
        Statements refining it in the hierarchy graph.

        Parameters
        ----------
        statements : list[indra.statements.Statement]
            A list of INDRA Statements whose belief scores are to
            be calculated. Each Statement object's belief attribute is updated
            by this function.
        """
        beliefs = self.get_hierarchy_probs(statements)
        for stmt in statements:
            sh = stmt.get_hash(self.matches_fun)
            stmt.belief = beliefs[sh]

    def get_hierarchy_probs(self, statements):
        # We only re-build the refinements graph if one wasn't provided
        # as an argument
        if self.refinements_graph is None:
            # Collect all the hashes and relevant statements, including
            # those that may not be in the given list but that are more
            # generic and in stmt.supported_by:
            stmts_by_hash = {}
            for stmt in statements:
                # The hash of *this* statement
                stmts_by_hash[stmt.get_hash(self.matches_fun)] = stmt
                for sb in stmt.supported_by:
                    stmts_by_hash[sb.get_hash(self.matches_fun)] = sb
            # Build the graph
            self.refinements_graph = \
                build_refinements_graph(stmts_by_hash=stmts_by_hash,
                                        matches_fun=self.matches_fun)
            assert_no_cycle(self.refinements_graph)
        logger.debug('Start belief calculation over refinements graph')
        # TODO: Change this to collect supporters data
        # structure (indexed by hash?) and score statements with the
        # accumulated evidence in one go.
        refiners_list = []
        # Note here that in collecting the list of statements that we pass
        # to get_refinement_probs, we build up a (potentially larger) list
        # of statements that includes the refined/supported_by statements of
        # the given statements
        graph_stmts = []
        for node in self.refinements_graph.nodes():
            # Get the statement
            stmt = self.refinements_graph.nodes[node]['stmt']
            graph_stmts.append(stmt)
            # Get the refiners/more specific stmts, if any
            refiners = list(networkx.descendants(self.refinements_graph, node))
            refiners_list.append(refiners)
        # Get dictionary mapping hashes to belief values
        beliefs_by_hash = self.get_refinement_probs(graph_stmts, refiners_list)
        logger.debug('Finished belief calculation over refinements graph')
        return beliefs_by_hash

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


def build_refinements_graph(stmts_by_hash, matches_fun=None):
    """Return a DiGraph based on matches hashes and Statement refinements.

    Parameters
    ----------
    stmts_by_hash : dict[int, indra.statements.Statement]
        A dict of statements keyed by their hashes.
    matches_fun : Optional[function]
        An optional function to calculate the matches key and hash of a
        given statement. Default: None

    Returns
    -------
    networkx.DiGraph
        A networkx graph whose nodes are statement hashes carrying a stmt
        attribute with the actual statement object. Edges point from
        less detailed to more detailed statements (i.e., from a statement
        to another statement that refines it).
    """
    logger.debug('Building refinements graph')
    g = networkx.DiGraph()
    for sh1, st1 in stmts_by_hash.items():
        g.add_node(sh1, stmt=st1)
        for st2 in st1.supported_by:
            sh2 = st2.get_hash(matches_fun=matches_fun)
            st2 = stmts_by_hash[sh2]
            g.add_node(sh2, stmt=st2)
            g.add_edge(sh2, sh1)
    logger.debug('Finished building refinements graph')
    return g


def extend_refinements_graph(g, stmt, less_specifics, matches_fun=None):
    """Extend refinements graph with a new statement and its refinements.

    Parameters
    ----------
    g : networkx.DiGraph
        A refinements graph to be extended.
    stmt : indra.statements.Statement
        The statement to be added to the refinements graph.
    less_specifics : list[int]
        A list of statement hashes of statements that are refined
        by this statement (i.e., are less specific versions of it).
    matches_fun : Optional[function]
        An optional function to calculate the matches key and hash of a
        given statement. Default: None
    """
    sh = stmt.get_hash(matches_fun=matches_fun)
    g.add_node(sh, stmt=stmt)
    for less_spec in less_specifics:
        g.add_edge(less_spec, sh)
    return g


def assert_no_cycle(g):
    """If the graph has cycles, throws AssertionError.

    This can be used to make sure that a refinements graph is a DAG.

    Parameters
    ----------
    g : networkx.DiGraph
        A refinements graph.
    """
    logger.debug('Looking for cycles in belief graph')
    try:
        cyc = networkx.algorithms.cycles.find_cycle(g)
    except networkx.exception.NetworkXNoCycle:
        return
    msg = 'Cycle found in hierarchy graph: %s' % cyc
    assert False, msg


def get_ranked_stmts(g):
    """Return a topological sort of statements from a graph."""
    logger.debug('Getting ranked statements')
    node_ranks = networkx.algorithms.dag.topological_sort(g)
    node_ranks = reversed(list(node_ranks))
    stmts = [g.nodes[n]['stmt'] for n in node_ranks]
    return stmts
