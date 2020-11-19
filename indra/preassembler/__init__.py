import time
import logging
import itertools
import functools
import collections
import networkx as nx
from indra.util import fast_deepcopy
from indra.statements import *
from indra.statements import stmt_type as indra_stmt_type

logger = logging.getLogger(__name__)


class Preassembler(object):
    """De-duplicates statements and arranges them in a specificity hierarchy.

    Parameters
    ----------
    ontology : :py:class:`indra.ontology.IndraOntology`
        An INDRA Ontology object.
    stmts : list of :py:class:`indra.statements.Statement` or None
        A set of statements to perform pre-assembly on. If None, statements
        should be added using the :py:meth:`add_statements` method.
    matches_fun : Optional[function]
        A functon which takes a Statement object as argument and
        returns a string key that is used for duplicate recognition. If
        supplied, it overrides the use of the built-in matches_key method of
        each Statement being assembled.
    refinement_fun : Optional[function]
        A function which takes two Statement objects and an ontology
        as an argument and returns True or False. If supplied, it overrides
        the built-in refinement_of method of each Statement being assembled.
    refinement_ns : Optional[set]
        A set of name spaces that should be considered for constructing
        refinements. If not provided, all name spaces are considered.
        Default: None

    Attributes
    ----------
    stmts : list of :py:class:`indra.statements.Statement`
        Starting set of statements for preassembly.
    unique_stmts : list of :py:class:`indra.statements.Statement`
        Statements resulting from combining duplicates.
    related_stmts : list of :py:class:`indra.statements.Statement`
        Top-level statements after building the refinement hierarchy.
    ontology : dict[:py:class:`indra.preassembler.ontology_graph.IndraOntology`]
        An INDRA Ontology object.
    """
    def __init__(self, ontology, stmts=None, matches_fun=None,
                 refinement_fun=None, refinement_ns=None):
        self.ontology = ontology
        if stmts:
            logger.debug("Deepcopying stmts in __init__")
            self.stmts = fast_deepcopy(stmts)
        else:
            self.stmts = []
        self.unique_stmts = None
        self.related_stmts = None
        self.matches_fun = matches_fun if matches_fun else \
            default_matches_fun
        self.refinement_fun = refinement_fun if refinement_fun else \
            default_refinement_fun
        self.refinement_ns = refinement_ns
        self._comparison_counter = 0

    def add_statements(self, stmts):
        """Add to the current list of statements.

        Parameters
        ----------
        stmts : list of :py:class:`indra.statements.Statement`
            Statements to add to the current list.
        """
        self.stmts += fast_deepcopy(stmts)

    def combine_duplicates(self):
        """Combine duplicates among `stmts` and save result in `unique_stmts`.

        A wrapper around the method :py:meth:`combine_duplicate_stmts`.
        """
        if self.unique_stmts is None:
            self.unique_stmts = self.combine_duplicate_stmts(self.stmts)
        return self.unique_stmts

    def _get_stmt_matching_groups(self, stmts):
        """Use the matches_fun method to get sets of matching statements."""
        # Remove exact duplicates using a set() call, then make copies:
        logger.debug('%d statements before removing object duplicates.' %
                     len(stmts))
        st = list(set(stmts))
        logger.debug('%d statements after removing object duplicates.' %
                     len(stmts))
        # Group statements according to whether they are matches (differing
        # only in their evidence).
        # Sort the statements in place by matches_key()
        st.sort(key=self.matches_fun)

        return itertools.groupby(st, key=self.matches_fun)

    def combine_duplicate_stmts(self, stmts):
        """Combine evidence from duplicate Statements.

        Statements are deemed to be duplicates if they have the same key
        returned by the `matches_key()` method of the Statement class. This
        generally means that statements must be identical in terms of their
        arguments and can differ only in their associated `Evidence` objects.

        This function keeps the first instance of each set of duplicate
        statements and merges the lists of Evidence from all of the other
        statements.

        Parameters
        ----------
        stmts : list of :py:class:`indra.statements.Statement`
            Set of statements to de-duplicate.

        Returns
        -------
        list of :py:class:`indra.statements.Statement`
            Unique statements with accumulated evidence across duplicates.

        Examples
        --------
        De-duplicate and combine evidence for two statements differing only
        in their evidence lists:

        >>> from indra.ontology.bio import bio_ontology
        >>> map2k1 = Agent('MAP2K1')
        >>> mapk1 = Agent('MAPK1')
        >>> stmt1 = Phosphorylation(map2k1, mapk1, 'T', '185',
        ... evidence=[Evidence(text='evidence 1')])
        >>> stmt2 = Phosphorylation(map2k1, mapk1, 'T', '185',
        ... evidence=[Evidence(text='evidence 2')])
        >>> pa = Preassembler(bio_ontology)
        >>> uniq_stmts = pa.combine_duplicate_stmts([stmt1, stmt2])
        >>> uniq_stmts
        [Phosphorylation(MAP2K1(), MAPK1(), T, 185)]
        >>> sorted([e.text for e in uniq_stmts[0].evidence])
        ['evidence 1', 'evidence 2']
        """
        # Helper function to get a list of evidence matches keys
        def _ev_keys(sts):
            ev_keys = []
            for stmt in sts:
                for ev in stmt.evidence:
                    ev_keys.append(ev.matches_key())
            return ev_keys
        # Iterate over groups of duplicate statements
        unique_stmts = []
        for _, duplicates in self._get_stmt_matching_groups(stmts):
            ev_keys = set()
            # Get the first statement and add the evidence of all subsequent
            # Statements to it
            duplicates = list(duplicates)
            start_ev_keys = _ev_keys(duplicates)
            for stmt_ix, stmt in enumerate(duplicates):
                if stmt_ix == 0:
                    new_stmt = stmt.make_generic_copy()
                if len(duplicates) == 1:
                    new_stmt.uuid = stmt.uuid
                raw_text = [None if ag is None else ag.db_refs.get('TEXT')
                            for ag in stmt.agent_list(deep_sorted=True)]
                raw_grounding = [None if ag is None else ag.db_refs
                                 for ag in stmt.agent_list(deep_sorted=True)]
                for ev in stmt.evidence:
                    ev_key = ev.matches_key() + str(raw_text) + \
                        str(raw_grounding)
                    if ev_key not in ev_keys:
                        # In case there are already agents annotations, we
                        # just add a new key for raw_text, otherwise create
                        # a new key
                        if 'agents' in ev.annotations:
                            ev.annotations['agents']['raw_text'] = raw_text
                            ev.annotations['agents']['raw_grounding'] = \
                                raw_grounding
                        else:
                            ev.annotations['agents'] = \
                                {'raw_text': raw_text,
                                 'raw_grounding': raw_grounding}
                        if 'prior_uuids' not in ev.annotations:
                            ev.annotations['prior_uuids'] = []
                        ev.annotations['prior_uuids'].append(stmt.uuid)
                        new_stmt.evidence.append(ev)
                        ev_keys.add(ev_key)
            end_ev_keys = _ev_keys([new_stmt])
            if len(end_ev_keys) != len(start_ev_keys):
                logger.debug('%d redundant evidences eliminated.' %
                             (len(start_ev_keys) - len(end_ev_keys)))
            # This should never be None or anything else
            assert isinstance(new_stmt, Statement)
            unique_stmts.append(new_stmt)
        return unique_stmts

    # Note that the kwargs here are just there for backwards compatibility
    # with old code that uses arguments related to multiprocessing.
    def combine_related(self, return_toplevel=True, filters=None, **kwargs):
        """Connect related statements based on their refinement relationships.

        This function takes as a starting point the unique statements (with
        duplicates removed) and returns a modified flat list of statements
        containing only those statements which do not represent a refinement of
        other existing statements. In other words, the more general versions of
        a given statement do not appear at the top level, but instead are
        listed in the `supports` field of the top-level statements.

        If :py:attr:`unique_stmts` has not been initialized with the
        de-duplicated statements, :py:meth:`combine_duplicates` is called
        internally.

        After this function is called the attribute :py:attr:`related_stmts` is
        set as a side-effect.

        The procedure for combining statements in this way involves a series
        of steps:

        1. The statements are grouped by type (e.g., Phosphorylation) and
           each type is iterated over independently.
        2. Each statement's agents are then aligned in a role-wise manner
           with the ontology being used, and all other statements which
           this statement can possibly refine are found.
        4. Each statement is then compared with the set of other statements
           identified earlier. If the statement represents a refinement of
           the other (as defined by the `refinement_of()` method implemented
           for the Statement), then the more refined statement is added
           to the `supports` field of the more general statement, and the
           more general statement is added to the `supported_by` field of
           the more refined statement.
        5. A new flat list of statements is created that contains only those
           statements that have no `supports` entries (statements containing
           such entries are not eliminated, because they will be retrievable
           from the `supported_by` fields of other statements). This list
           is returned to the caller.

        .. note:: Subfamily relationships must be consistent across arguments

            For now, we require that merges can only occur if the *isa*
            relationships are all in the *same direction for all the agents* in
            a Statement. For example, the two statement groups: `RAF_family ->
            MEK1` and `BRAF -> MEK_family` would not be merged, since BRAF
            *isa* RAF_family, but MEK_family is not a MEK1. In the future this
            restriction could be revisited.

        Parameters
        ----------
        return_toplevel : Optional[bool]
            If True only the top level statements are returned.
            If False, all statements are returned. Default: True
        filters : Optional[list[function]]
            A list of function handles that define filter functions on
            possible statement refinements. Each function takes
            a stmts_by_hash dict and a stmts_to_compare dict as its input and
            returns a dict of possible refinements where the keys are
            statement hashes and the values are sets of statement hashes that
            the key statement possibly refines. If not provided, a built-in
            ontology-based pre-filter is applied. Note, that if a list of filter
            functions is provided, the built-in ontology-based pre-filter is not
            automatically appended to the list of filters. In this case,
            consider adding the `ontology_refinement_filter` function from this
            module to the filters list.

        Returns
        -------
        list of :py:class:`indra.statement.Statement`
            The returned list contains Statements representing the more
            concrete/refined versions of the Statements involving particular
            entities. The attribute :py:attr:`related_stmts` is also set to
            this list. However, if return_toplevel is False then all
            statements are returned, irrespective of level of specificity.
            In this case the relationships between statements can
            be accessed via the supports/supported_by attributes.

        Examples
        --------
        A more general statement with no information about a Phosphorylation
        site is identified as supporting a more specific statement:

        >>> from indra.ontology.bio import bio_ontology
        >>> braf = Agent('BRAF')
        >>> map2k1 = Agent('MAP2K1')
        >>> st1 = Phosphorylation(braf, map2k1)
        >>> st2 = Phosphorylation(braf, map2k1, residue='S')
        >>> pa = Preassembler(bio_ontology, [st1, st2])
        >>> combined_stmts = pa.combine_related() # doctest:+ELLIPSIS
        >>> combined_stmts
        [Phosphorylation(BRAF(), MAP2K1(), S)]
        >>> combined_stmts[0].supported_by
        [Phosphorylation(BRAF(), MAP2K1())]
        >>> combined_stmts[0].supported_by[0].supports
        [Phosphorylation(BRAF(), MAP2K1(), S)]
        """
        if self.related_stmts is not None:
            if return_toplevel:
                return self.related_stmts
            else:
                assert self.unique_stmts is not None
                return self.unique_stmts

        # Call combine_duplicates, which lazily initializes self.unique_stmts
        unique_stmts = self.combine_duplicates()

        # Generate the index map, linking related statements.
        idx_map = self._generate_id_maps(unique_stmts,
                                         filters=filters)

        # Now iterate over all indices and set supports/supported by
        for ix1, ix2 in idx_map:
            unique_stmts[ix1].supported_by.append(unique_stmts[ix2])
            unique_stmts[ix2].supports.append(unique_stmts[ix1])
        # Get the top level statements
        self.related_stmts = [st for st in unique_stmts if not st.supports]
        logger.debug('%d top level' % len(self.related_stmts))
        if return_toplevel:
            return self.related_stmts
        else:
            return unique_stmts

    # Note that the kwargs here are just there for backwards compatibility
    # with old code that uses arguments related to multiprocessing.
    def _generate_id_maps(self, unique_stmts, split_idx=None,
                          filters=None, **kwargs):
        """Return pairs of statement indices representing refinement relations.

        Parameters
        ----------
        unique_stmts : list[indra.statements.Statement]
            A list of de-duplicated INDRA Statements.
        split_idx : Optional[int]
            An index at which the flat list of unique statements should be split
            and compared for refinements only across the two groups, not
            within each group.
        filters : Optional[list[function]]
            A list of function handles that define filter functions on
            possible statement refinements. Each function takes
            a stmts_by_hash dict as its input and returns a dict
            of possible refinements where the keys are statement hashes
            and the values are sets of statement hashes that the
            key statement possibly refines.

        Returns
        -------
        list[tuple]
            A list of tuples where the first element of each tuple is
            the linear index of a statement in the unique stmts list
            which refines the statement whose index is the second
            element of the tuple.
        """
        ts = time.time()
        # Make a list of Statement types
        stmt_to_idx = {stmt.get_hash(matches_fun=self.matches_fun): idx
                       for idx, stmt in enumerate(unique_stmts)}

        # Statements keyed by their hashes
        stmts_by_hash = {stmt.get_hash(matches_fun=self.matches_fun):
                         stmt for stmt in unique_stmts}
        stmts_to_compare = None
        # Here we apply any additional filters to cut down the number of
        # potential comparisons before actually making comparisons
        if filters:
            # We apply filter functions sequentially
            for filter_fun in filters:
                logger.debug('Applying filter %s' % filter_fun.__name__)
                stmts_to_compare = \
                    filter_fun(stmts_by_hash, stmts_to_compare)
                total_comparisons = sum(len(v)
                                        for v in stmts_to_compare.values())
                logger.debug('Total comparisons after filter %s: %d' %
                             (filter_fun.__name__, total_comparisons))
        else:
            stmts_to_compare = \
                ontology_refinement_filter(stmts_by_hash=stmts_by_hash,
                                           stmts_to_compare=stmts_to_compare,
                                           ontology=self.ontology)
            total_comparisons = sum(len(v) for v in stmts_to_compare.values())
            logger.debug('Total comparisons after ontology filter: %d' %
                         total_comparisons)

        te = time.time()
        logger.info('Applied all refinement pre-filters in %.2fs' % (te-ts))
        logger.info('Total comparisons: %d' % total_comparisons)

        # Here we handle split_idx to allow finding refinements between
        # to distinct groups of statements (identified by an index at which we
        # split the unique_statements list) rather than globally across
        # all unique statements.
        if split_idx:
            # This dict maps statement hashes to a bool value based on which
            # of the two groups the statement belongs to.
            hash_to_split_group = {sh: (idx <= split_idx) for sh, idx
                                   in stmt_to_idx.items()}
        else:
            hash_to_split_group = None

        # We can now do the actual comparisons and return pairs of confirmed
        # refinements in a list.
        maps = \
            self.confirm_possible_refinements(stmts_by_hash,
                                              stmts_to_compare,
                                              split_groups=hash_to_split_group)

        idx_maps = [(stmt_to_idx[refinement], stmt_to_idx[refined])
                    for refinement, refined in maps]
        return idx_maps

    def confirm_possible_refinements(self, stmts_by_hash, stmts_to_compare,
                                     split_groups=None):
        """Return confirmed pairs of statement refinement relationships.

        Parameters
        ----------
        stmts_by_hash : dict
            A dict whose keys are statement hashes that point to the
            (deduplicated) statement with that hash as a value.
        stmts_to_compare : dict
            A dict whose keys are statement hashes and values are sets of
            statement hashes that the statement with the given hash can
            possibly refine.
        split_groups : dict
            A dict whose keys are statement hashes and values represent
            one of two groups that the statement is in. Statement in the
            same group aren't compared, only statements in different
            groups are. This can be used to do "bipartite" refinement
            checking across a set of statements.

        Returns
        -------
        list of tuple
            A list of tuple where the first element of each tuple is the
            hash of a statement which refines that statement whose hash
            is the second element of the tuple.
        """
        maps = []
        # We again iterate over statements
        ts = time.time()
        # Given the possible refinements in stmts_to_compare, we confirm each
        for stmt_hash, possible_refined_hashes in stmts_to_compare.items():
            # We use the previously constructed set of statements that this one
            # can possibly refine
            for possible_refined_hash in possible_refined_hashes:
                # We handle split groups here to only check refinements between
                # statements that are in different groups to compare
                if not split_groups or split_groups[stmt_hash] != \
                        split_groups[possible_refined_hash]:
                    # And then do the actual comparison. Here we use
                    # entities_refined=True which means that we assert that
                    # the entities, in each role, are already confirmed to
                    # be "compatible" for refinement, and therefore, we
                    # don't need to again confirm this (i.e., call "isa") in
                    # the refinement_of function.
                    self._comparison_counter += 1
                    ref = self.refinement_fun(
                        stmts_by_hash[stmt_hash],
                        stmts_by_hash[possible_refined_hash],
                        ontology=self.ontology,
                        # NOTE: here we assume that the entities at this point
                        # are definitely refined due to the use of an
                        # ontology-based pre-filter. If this is not the case
                        # for some reason then it is the responsibility of the
                        # user-supplied self.refinement_fun to disregard the
                        # entities_refined argument.
                        entities_refined=True)
                    if ref:
                        maps.append((stmt_hash, possible_refined_hash))
        te = time.time()
        logger.debug('Confirmed %d refinements in %.2fs' % (len(maps), te-ts))
        return maps

    def find_contradicts(self):
        """Return pairs of contradicting Statements.

        Returns
        -------
        contradicts : list(tuple(Statement, Statement))
            A list of Statement pairs that are contradicting.
        """
        # Make a dict of Statement by type
        stmts_by_type = collections.defaultdict(list)
        for stmt in self.stmts:
            stmts_by_type[indra_stmt_type(stmt)].append(stmt)
        stmts_by_type = dict(stmts_by_type)

        # Handle Statements with polarity first
        pos_stmts = AddModification.__subclasses__()
        neg_stmts = [modclass_to_inverse[c] for c in pos_stmts]

        pos_stmts += [Activation, IncreaseAmount]
        neg_stmts += [Inhibition, DecreaseAmount]

        contradicts = []
        # Handle statements with polarity first
        # TODO: we could probably do some optimization here
        # to not have to check statements combinatorially
        for pst, nst in zip(pos_stmts, neg_stmts):
            poss = stmts_by_type.get(pst, [])
            negs = stmts_by_type.get(nst, [])

            for ps, ns in itertools.product(poss, negs):
                if ps.contradicts(ns, self.ontology):
                    contradicts.append((ps, ns))

        # Handle neutral Statements next
        neu_stmts = [Influence, ActiveForm]
        for stt in neu_stmts:
            stmts = stmts_by_type.get(stt, [])
            for st1, st2 in itertools.combinations(stmts, 2):
                if st1.contradicts(st2, self.ontology):
                    contradicts.append((st1, st2))

        return contradicts

    def _normalize_relations(self, ns, rank_key, rel_fun, flip_polarity):
        # Find related entries, sort them, and return the first one which is
        # the one that will be normalized to
        def _replace_grounding(ns, entry, rank_key, rel_fun):
            rel_ents = rel_fun(ns, entry)
            if rel_ents:
                rel_ents = [(ns, e.split('#')[1] if '#' in e else e)
                            for ns, e in rel_ents]
                sorted_entries = sorted([(ns, entry)] + rel_ents,
                                        key=rank_key)
                _, chosen = sorted_entries[0]
                return chosen, chosen != entry
            else:
                return entry, False

        # If no custom rank_key was provided we use the original value to
        # sort by
        if rank_key is None:
            def polarity_rank_key(args):
                ns, entry = args
                pol = self.ontology.get_polarity(ns, entry)
                # Here we flip polarities to rank positive polarity before
                # negative
                pol_rank = -1 if pol is None else -pol
                return pol_rank, entry
            rank_key = polarity_rank_key
        # We now go agent by agent to normalize grounding
        for stmt in self.stmts:
            for agent_idx, agent in enumerate(stmt.agent_list()):
                # If the relevant namespace is an entry
                if agent is not None and ns in agent.db_refs:
                    grounding = agent.db_refs[ns]
                    # If we have a list, we iterate over it and normalize
                    # each entry separately
                    if isinstance(grounding, list):
                        new_grounding = []
                        for idx, (entry, score) in enumerate(grounding):
                            chosen, changed = _replace_grounding(ns, entry,
                                                                 rank_key,
                                                                 rel_fun)
                            new_grounding.append((chosen, score))
                            # If the top grounding was changed and we need
                            # to flip polarity then the Statement's polarity
                            # is flipped
                            if idx == 0 and changed and flip_polarity:
                                stmt.flip_polarity(agent_idx=agent_idx)
                        agent.db_refs[ns] = new_grounding
                    # If there's only one grounding then we just normalize
                    # that one
                    else:
                        chosen, changed = _replace_grounding(ns, grounding,
                                                             rank_key, rel_fun)
                        agent.db_refs[ns] = chosen
                        if changed and flip_polarity:
                            stmt.flip_polarity(agent_idx=agent_idx)

    def normalize_equivalences(self, ns, rank_key=None):
        """Normalize to one of a set of equivalent concepts across statements.

        This function changes Statements in place without returning a value.

        Parameters
        ----------
        ns : str
            The db_refs namespace for which the equivalence relation should
            be applied.
        rank_key : Optional[function]
            A function handle which assigns a sort key to each entry in the
            given namespace to allow prioritizing in a controlled way which
            concept is normalized to.
        """
        rel_fun = functools.partial(self.ontology.child_rel,
                                    rel_types={'is_equal'})
        self._normalize_relations(ns, rank_key, rel_fun, False)

    def normalize_opposites(self, ns, rank_key=None):
        """Normalize to one of a pair of opposite concepts across statements.

        This function changes Statements in place without returning a value.

        Parameters
        ----------
        ns : str
            The db_refs namespace for which the opposite relation should
            be applied.
        rank_key : Optional[function]
            A function handle which assigns a sort key to each entry in the
            given namespace to allow prioritizing in a controlled way which
            concept is normalized to.
        """
        rel_fun = functools.partial(self.ontology.child_rel,
                                    rel_types={'is_opposite'})
        self._normalize_relations(ns, rank_key, rel_fun, True)


def render_stmt_graph(statements, reduce=True, english=False, rankdir=None,
                      agent_style=None):
    """Render the statement hierarchy as a pygraphviz graph.

    Parameters
    ----------
    statements : list of :py:class:`indra.statements.Statement`
        A list of top-level statements with associated supporting statements
        resulting from building a statement hierarchy with
        :py:meth:`combine_related`.
    reduce : bool
        Whether to perform a transitive reduction of the edges in the graph.
        Default is True.
    english : bool
        If True, the statements in the graph are represented by their
        English-assembled equivalent; otherwise they are represented as
        text-formatted Statements.
    rankdir : str or None
        Argument to pass through to the  pygraphviz `AGraph` constructor
        specifying graph layout direction. In particular, a value of 'LR'
        specifies a left-to-right direction. If None, the pygraphviz default
        is used.
    agent_style : dict or None
        Dict of attributes specifying the visual properties of nodes. If None,
        the following default attributes are used::

            agent_style = {'color': 'lightgray', 'style': 'filled',
                           'fontname': 'arial'}

    Returns
    -------
    pygraphviz.AGraph
        Pygraphviz graph with nodes representing statements and edges pointing
        from supported statements to supported_by statements.

    Examples
    --------
    Pattern for getting statements and rendering as a Graphviz graph:

    >>> from indra.ontology.bio import bio_ontology
    >>> braf = Agent('BRAF')
    >>> map2k1 = Agent('MAP2K1')
    >>> st1 = Phosphorylation(braf, map2k1)
    >>> st2 = Phosphorylation(braf, map2k1, residue='S')
    >>> pa = Preassembler(bio_ontology, [st1, st2])
    >>> pa.combine_related() # doctest:+ELLIPSIS
    [Phosphorylation(BRAF(), MAP2K1(), S)]
    >>> graph = render_stmt_graph(pa.related_stmts)
    >>> graph.write('example_graph.dot') # To make the DOT file
    >>> graph.draw('example_graph.png', prog='dot') # To make an image

    Resulting graph:

    .. image:: /images/example_graph.png
        :align: center
        :alt: Example statement graph rendered by Graphviz

    """
    import pygraphviz as pgv
    from indra.assemblers.english import EnglishAssembler
    # Set the default agent formatting properties
    if agent_style is None:
        agent_style = {'color': 'lightgray', 'style': 'filled',
                       'fontname': 'arial'}
    # Sets to store all of the nodes and edges as we recursively process all
    # of the statements
    nodes = set([])
    edges = set([])
    stmt_dict = {}

    # Recursive function for processing all statements
    def process_stmt(stmt):
        nodes.add(str(stmt.matches_key()))
        stmt_dict[str(stmt.matches_key())] = stmt
        for sby_ix, sby_stmt in enumerate(stmt.supported_by):
            edges.add((str(stmt.matches_key()), str(sby_stmt.matches_key())))
            process_stmt(sby_stmt)

    # Process all of the top-level statements, getting the supporting statements
    # recursively
    for stmt in statements:
        process_stmt(stmt)
    # Create a networkx graph from the nodes
    nx_graph = nx.DiGraph()
    nx_graph.add_edges_from(edges)
    # Perform transitive reduction if desired
    if reduce:
        nx_graph = nx.algorithms.dag.transitive_reduction(nx_graph)
    # Create a pygraphviz graph from the nx graph
    try:
        pgv_graph = pgv.AGraph(name='statements', directed=True,
                               rankdir=rankdir)
    except NameError:
        logger.error('Cannot generate graph because '
                     'pygraphviz could not be imported.')
        return None
    for node in nx_graph.nodes():
        stmt = stmt_dict[node]
        if english:
            ea = EnglishAssembler([stmt])
            stmt_str = ea.make_model()
        else:
            stmt_str = str(stmt)
        pgv_graph.add_node(node,
                           label='%s (%d)' % (stmt_str, len(stmt.evidence)),
                           **agent_style)
    pgv_graph.add_edges_from(nx_graph.edges())
    return pgv_graph


def flatten_stmts(stmts):
    """Return the full set of unique stms in a pre-assembled stmt graph.

    The flattened list of statements returned by this function can be
    compared to the original set of unique statements to make sure no
    statements have been lost during the preassembly process.

    Parameters
    ----------
    stmts : list of :py:class:`indra.statements.Statement`
        A list of top-level statements with associated supporting statements
        resulting from building a statement hierarchy with
        :py:meth:`combine_related`.

    Returns
    -------
    stmts : list of :py:class:`indra.statements.Statement`
        List of all statements contained in the hierarchical statement graph.

    Examples
    --------
    Calling :py:meth:`combine_related` on two statements results in one
    top-level statement; calling :py:func:`flatten_stmts` recovers both:

    >>> from indra.ontology.bio import bio_ontology
    >>> braf = Agent('BRAF')
    >>> map2k1 = Agent('MAP2K1')
    >>> st1 = Phosphorylation(braf, map2k1)
    >>> st2 = Phosphorylation(braf, map2k1, residue='S')
    >>> pa = Preassembler(bio_ontology, [st1, st2])
    >>> pa.combine_related() # doctest:+ELLIPSIS
    [Phosphorylation(BRAF(), MAP2K1(), S)]
    >>> flattened = flatten_stmts(pa.related_stmts)
    >>> flattened.sort(key=lambda x: x.matches_key())
    >>> flattened
    [Phosphorylation(BRAF(), MAP2K1()), Phosphorylation(BRAF(), MAP2K1(), S)]
    """
    total_stmts = set(stmts)
    for stmt in stmts:
        if stmt.supported_by:
            children = flatten_stmts(stmt.supported_by)
            total_stmts = total_stmts.union(children)
    return list(total_stmts)


def flatten_evidence(stmts, collect_from=None):
    """Add evidence from *supporting* stmts to evidence for *supported* stmts.

    Parameters
    ----------
    stmts : list of :py:class:`indra.statements.Statement`
        A list of top-level statements with associated supporting statements
        resulting from building a statement hierarchy with
        :py:meth:`combine_related`.
    collect_from : str in ('supports', 'supported_by')
        String indicating whether to collect and flatten evidence from the
        `supports` attribute of each statement or the `supported_by` attribute.
        If not set, defaults to 'supported_by'.

    Returns
    -------
    stmts : list of :py:class:`indra.statements.Statement`
        Statement hierarchy identical to the one passed, but with the
        evidence lists for each statement now containing all of the evidence
        associated with the statements they are supported by.

    Examples
    --------
    Flattening evidence adds the two pieces of evidence from the supporting
    statement to the evidence list of the top-level statement:

    >>> from indra.ontology.bio import bio_ontology
    >>> braf = Agent('BRAF')
    >>> map2k1 = Agent('MAP2K1')
    >>> st1 = Phosphorylation(braf, map2k1,
    ... evidence=[Evidence(text='foo'), Evidence(text='bar')])
    >>> st2 = Phosphorylation(braf, map2k1, residue='S',
    ... evidence=[Evidence(text='baz'), Evidence(text='bak')])
    >>> pa = Preassembler(bio_ontology, [st1, st2])
    >>> pa.combine_related() # doctest:+ELLIPSIS
    [Phosphorylation(BRAF(), MAP2K1(), S)]
    >>> [e.text for e in pa.related_stmts[0].evidence]
    ['baz', 'bak']
    >>> flattened = flatten_evidence(pa.related_stmts)
    >>> sorted([e.text for e in flattened[0].evidence])
    ['bak', 'bar', 'baz', 'foo']
    """
    if collect_from is None:
        collect_from = 'supported_by'
    if collect_from not in ('supports', 'supported_by'):
        raise ValueError('collect_from must be one of "supports", '
                         '"supported_by"')
    logger.info('Flattening evidence based on %s' % collect_from)
    # Copy all of the statements--these will be the ones where we update
    # the evidence lists
    stmts = fast_deepcopy(stmts)
    for stmt in stmts:
        # We get the original evidence keys here so we can differentiate them
        # from ones added during flattening.
        orig_ev_keys = [ev.matches_key() for ev in stmt.evidence]
        # We now do the flattening
        total_evidence = _flatten_evidence_for_stmt(stmt, collect_from)
        # Here we add annotations for each evidence in the list,
        # depending on whether it's an original direct evidence or one that
        # was added during flattening
        new_evidence = []
        for ev in total_evidence:
            ev_key = ev.matches_key()
            if ev_key in orig_ev_keys:
                ev.annotations['support_type'] = 'direct'
                new_evidence.append(ev)
            else:
                ev_copy = fast_deepcopy(ev)
                ev_copy.annotations['support_type'] = collect_from
                new_evidence.append(ev_copy)
        # Now set the new evidence list as the copied statement's evidence
        stmt.evidence = new_evidence
    return stmts


def _flatten_evidence_for_stmt(stmt, collect_from):
    supp_stmts = (stmt.supports if collect_from == 'supports'
                  else stmt.supported_by)
    total_evidence = set(stmt.evidence)
    for supp_stmt in supp_stmts:
        child_evidence = _flatten_evidence_for_stmt(supp_stmt, collect_from)
        total_evidence = total_evidence.union(child_evidence)
    return list(total_evidence)


def default_refinement_fun(st1, st2, ontology, entities_refined):
    return st1.refinement_of(st2, ontology, entities_refined)


def default_matches_fun(st):
    return st.matches_key()


# TODO: we could make the agent key function parameterizable with the
# preassembler to allow custom agent mappings to the ontology.
def get_agent_key(agent):
    """Return a key for an Agent for use in refinement finding.

    Parameters
    ----------
    agent : indra.statements.Agent or None
         An INDRA Agent whose key should be returned.

    Returns
    -------
    tuple or None
        The key that maps the given agent to the ontology, with special
        handling for ungrounded and None Agents.
    """
    if isinstance(agent, Event):
        agent = agent.concept
    if agent is None:
        agent_key = None
    else:
        agent_key = agent.get_grounding()
        if not agent_key[0]:
            agent_key = ('NAME', agent.name)
    return agent_key


def get_relevant_keys(agent_key, all_keys_for_role, ontology):
    """Return relevant agent keys for an agent key for refinement finding.

    Parameters
    ----------
    agent_key : tuple or None
        An agent key of interest.
    all_keys_for_role : set
        The set of all agent keys in a given statement corpus with a
        role matching that of the given agent_key.
    ontology : indra.ontology.IndraOntology
        An IndraOntology instance with respect to which relevant other
        agent keys are found for the purposes of refinement.

    Returns
    -------
    set
        The set of relevant agent keys which this given agent key can
        possibly refine.
    """
    relevant_keys = {None, agent_key}
    if agent_key is not None:
        relevant_keys |= set(ontology.get_parents(*agent_key))
    relevant_keys &= all_keys_for_role
    return relevant_keys


def ontology_refinement_filter(stmts_by_hash, stmts_to_compare, ontology):
    """Return possible refinement relationships based on an ontology.

    Parameters
    ----------
    stmts_by_hash : dict
        A dict whose keys are statement hashes that point to the
        (deduplicated) statement with that hash as a value.
    stmts_to_compare : dict or None
        A dict of existing statements to compare that will be further
        filtered down in this function and then returned.
    ontology : indra.ontology.IndraOntology
        An IndraOntology instance iwth respect to which this
        filter is applied.

    Returns
    -------
    dict
        A dict whose keys are statement hashes and values are sets
        of statement hashes that can potentially be refined by the
        statement identified by the key.
    """
    ts = time.time()
    stmts_by_type = collections.defaultdict(set)
    for stmt_hash, stmt in stmts_by_hash.items():
        stmts_by_type[indra_stmt_type(stmt)].add(stmt_hash)
    stmts_by_type = dict(stmts_by_type)

    for stmt_type, stmt_hashes in stmts_by_type.items():
        logger.info('Finding ontology-based refinements for %d %s statements'
                    % (len(stmts_by_hash), stmt_type.__name__))
        stmts_by_hash_this_type = {
            stmt_hash: stmts_by_hash[stmt_hash]
            for stmt_hash in stmt_hashes
        }
        stmts_to_compare = \
            ontology_refinement_filter_by_stmt_type(stmts_by_hash_this_type,
                                                    stmts_to_compare,
                                                    ontology)
    te = time.time()
    logger.debug('Identified ontology-based possible refinements in %.2fs'
                 % (te-ts))
    # Make an empty dict to make sure we don't return a None
    if stmts_to_compare is None:
        stmts_to_compare = {}
    return stmts_to_compare


def ontology_refinement_filter_by_stmt_type(stmts_by_hash, stmts_to_compare,
                                            ontology):
    """Return possible refinement relationships based on an ontology.

    Importantly, here we assume that all statements in stmts_by_hash
    are of a single type.

    Parameters
    ----------
    stmts_by_hash : dict
        A dict whose keys are statement hashes that point to the
        (deduplicated) statement with that hash as a value.
    stmts_to_compare : dict or None
        A dict of existing statements to compare that will be further
        filtered down in this function and then returned.
    ontology : indra.ontology.IndraOntology
        An IndraOntology instance iwth respect to which this
        filter is applied.

    Returns
    -------
    list of tuple
        A list of tuples where the first element of each tuple is the
        hash of a statement which refines that statement whose hash
        is the second element of the tuple.
    """
    # Step 1. initialize data structures
    roles = stmts_by_hash[next(iter(stmts_by_hash))]._agent_order
    # Mapping agent keys to statement hashes
    agent_key_to_hash = {}
    # Mapping statement hashes to agent keys
    hash_to_agent_key = {}
    # All agent keys for a given agent role
    all_keys_by_role = {}
    for role in roles:
        agent_key_to_hash[role] = collections.defaultdict(set)
        hash_to_agent_key[role] = collections.defaultdict(set)

    # Step 2. Fill up the initial data structures in preparation
    # for identifying potential refinements
    for sh, stmt in stmts_by_hash.items():
        for role in roles:
            agents = getattr(stmt, role)
            # Handle a special case here where a list=like agent
            # role can be empty, here we will consider anything else
            # to be a refinement, hence add a None key
            if isinstance(agents, list) and not agents:
                agent_keys = {None}
            # Generally, we take all the agent keys for a single or
            # list-like agent role.
            else:
                agent_keys = {get_agent_key(agent) for agent in
                              (agents if isinstance(agents, list)
                               else [agents])}
            for agent_key in agent_keys:
                agent_key_to_hash[role][agent_key].add(sh)
                hash_to_agent_key[role][sh].add(agent_key)

    agent_key_to_hash = dict(agent_key_to_hash)
    hash_to_agent_key = dict(hash_to_agent_key)

    for role in roles:
        all_keys_by_role[role] = set(agent_key_to_hash[role].keys())

    # Step 3. Identify all the pairs of statements which can be in a
    # refinement relationship
    first_filter = True if stmts_to_compare is None else False
    if first_filter:
        stmts_to_compare = {}
    # We iterate over each statement and find all other statements that it
    # can potentially refine
    ts = time.time()
    for sh, stmt in stmts_by_hash.items():
        relevants = None
        # We now iterate over all the agent roles in the given statement
        # type
        for role, hash_to_agent_key_for_role in hash_to_agent_key.items():
            # We get all the agent keys in all other statements that the
            # agent
            # in this role in this statement can be a refinement.
            for agent_key in hash_to_agent_key_for_role[sh]:
                relevant_keys = get_relevant_keys(
                    agent_key,
                    all_keys_by_role[role],
                    ontology)
                # We now get the actual statement hashes that these other
                # potentially refined agent keys appear in in the given role
                role_relevant_stmt_hashes = set.union(
                    *[agent_key_to_hash[role][rel]
                      for rel in relevant_keys]) - {sh}
                # In the first iteration, we initialize the set with the
                # relevant statement hashes
                if relevants is None:
                    relevants = role_relevant_stmt_hashes
                # In subsequent iterations, we take the intersection of
                # the relevant sets per role
                else:
                    relevants &= role_relevant_stmt_hashes
        # These hashes are now the ones that this statement needs
        # to be compared against. Importantly, the relationship is in
        # a well-defined direction so we don't need to test both ways.
        if first_filter:
            stmts_to_compare[sh] = relevants
        else:
            stmts_to_compare[sh] = \
                stmts_to_compare.get(sh, set()) & relevants
    return stmts_to_compare


def bio_ontology_refinement_filter(stmts_by_hash, stmts_to_compare):
    """An ontology refinement filter that works with the INDRA BioOntology."""
    from indra.ontology.bio import bio_ontology
    return ontology_refinement_filter(stmts_by_hash, stmts_to_compare,
                                      ontology=bio_ontology)
