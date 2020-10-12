import sys
import time
import logging
import itertools
import functools
import collections
import networkx as nx
import multiprocessing as mp
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
                if stmt_ix is 0:
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

    def _generate_id_maps(self, unique_stmts, poolsize=None,
                          size_cutoff=100, split_idx=None):
        """Connect statements using their refinement relationship."""
        # FIXME: we probably want to use the custom matches key function here
        # Make a list of Statement types
        stmts_by_type = collections.defaultdict(list)
        stmt_to_idx = {stmt.get_hash(): idx
                       for idx, stmt in enumerate(unique_stmts)}
        for idx, stmt in enumerate(unique_stmts):
            stmts_by_type[indra_stmt_type(stmt)].append(stmt)

        maps = []
        import tqdm
        for stmts in tqdm.tqdm(stmts_by_type.values()):
            maps += self._generate_hash_maps_by_stmt_type(
                stmts, stmts[0]._agent_order)
        idx_maps = [(stmt_to_idx[refinement], stmt_to_idx[refined])
                    for refinement, refined in maps]
        return idx_maps

    def _generate_hash_maps_by_stmt_type(self, stmts, roles):
        # Step 1. initialize data structures
        # Statements keyed by their hashes
        stmts_by_hash = {stmt.get_hash(): stmt for stmt in stmts}
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
        # for ideentifying potential refinements
        for sh, stmt in stmts_by_hash.items():
            for role in roles:
                agents = getattr(stmt, role)
                agent_keys = get_agent_keys(agents)
                for agent_key in agent_keys:
                    agent_key_to_hash[role][agent_key].add(sh)
                hash_to_agent_key[role][sh].add(agent_key)

        for role in roles:
            all_keys_by_role[role] = set(agent_key_to_hash[role].keys())

        # Step 3. Identify all the pairs of statements between which can be
        # in a refinement relationship
        stmts_to_compare = {}
        # We iterate over each statement and find all other statements that it
        # can potentially refine
        for sh, stmt in stmts_by_hash.items():
            relevants = None
            # We now iterate over all the agent roles in the given statement
            # type
            for role, hash_to_agent_key_for_role in hash_to_agent_key.items():
                # We get all the agent keys in all other statements that the
                # agent
                # in this role in this statement can be a refinement.
                relevant_keys = get_relevant_keys(
                    hash_to_agent_key_for_role[sh],
                    all_keys_by_role[role],
                    self.ontology)
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
            stmts_to_compare[sh] = relevants

        total_comparisons = sum(len(v) for v in stmts_to_compare.values())
        logger.info('Total comparisons: %d' % total_comparisons)

        # Step 4. We can now do the actual comparisons and save pairs of
        # confirmed refinements in a list.
        maps = []
        # We again iterate over statements
        for stmt_hash, possible_refineds in stmts_to_compare.items():
            # We use the previously constructed set of statements that this one
            # can possibly refine
            for possible_refined in possible_refineds:
                # FIXME: custom functions?
                # And then do the actual comparison. Here we use
                # entities_refined=True which means that we assert that
                # the entities, in each role, are already confirmed to
                # be "compatible" for refinement, and therefore, we don't need
                # to again confirm this (i.e., call "isa") in the refinement_of
                # function.
                if stmts_by_hash[stmt_hash].refinement_of(
                        stmts_by_hash[possible_refined], self.ontology,
                        entities_refined=True):
                    maps.append((stmt_hash, possible_refined))
        return maps

    def combine_related(self, return_toplevel=True, poolsize=None,
                        size_cutoff=100):
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

        On multi-core machines, the algorithm can be parallelized by setting
        the poolsize argument to the desired number of worker processes.
        This feature is only available in Python > 3.4.

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
        poolsize : Optional[int]
            The number of worker processes to use to parallelize the
            comparisons performed by the function. If None (default), no
            parallelization is performed. NOTE: Parallelization is only
            available on Python 3.4 and above.
        size_cutoff : Optional[int]
            Groups with size_cutoff or more statements are sent to worker
            processes, while smaller groups are compared in the parent process.
            Default value is 100. Not relevant when parallelization is not
            used.

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
        idx_map = self._generate_id_maps(unique_stmts, poolsize, size_cutoff)

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

    def find_contradicts(self):
        """Return pairs of contradicting Statements.

        Returns
        -------
        contradicts : list(tuple(Statement, Statement))
            A list of Statement pairs that are contradicting.
        """
        # Make a dict of Statement by type
        stmts_by_type = collections.defaultdict(lambda: [])
        for idx, stmt in enumerate(self.stmts):
            stmts_by_type[indra_stmt_type(stmt)].append((idx, stmt))

        # Handle Statements with polarity first
        pos_stmts = AddModification.__subclasses__()
        neg_stmts = [modclass_to_inverse[c] for c in pos_stmts]

        pos_stmts += [Activation, IncreaseAmount]
        neg_stmts += [Inhibition, DecreaseAmount]

        contradicts = []
        for pst, nst in zip(pos_stmts, neg_stmts):
            poss = stmts_by_type.get(pst, [])
            negs = stmts_by_type.get(nst, [])

            pos_stmt_by_group = self._get_stmt_by_group(pst, poss,
                                                        self.ontology)
            neg_stmt_by_group = self._get_stmt_by_group(nst, negs,
                                                        self.ontology)
            for key, pg in pos_stmt_by_group.items():
                ng = neg_stmt_by_group.get(key, [])
                for (_, st1), (_, st2) in itertools.product(pg, ng):
                    if st1.contradicts(st2, self.ontology):
                        contradicts.append((st1, st2))

        # Handle neutral Statements next
        neu_stmts = [Influence, ActiveForm]
        for stt in neu_stmts:
            stmts = stmts_by_type.get(stt, [])
            for (_, st1), (_, st2) in itertools.combinations(stmts, 2):
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
    stmts : list of :py:class:`indra.statements.Statement`
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
    rank_dir : str or None
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


def default_refinement_fun(st1, st2, ontology):
    return st1.refinement_of(st2, ontology)


def default_matches_fun(st):
    return st.matches_key()


def get_agent_keys(agents):
    # To
    if not isinstance(agents, list):
        agents = {agents}
    keys = set()
    for agent in agents:
        if isinstance(agent, Event):
            agent = agent.concept
        if agent is None:
            agent_key = None
        else:
            agent_key = agent.get_grounding()
            if not agent_key[0]:
                agent_key = ('NAME', agent.name)
        keys.add(agent_key)
    return keys


def get_relevant_keys(agent_keys, all_keys_for_role, ontology):
    relevant_keys = {None}
    for agent_key in agent_keys:
        relevant_keys |= {agent_key}
        if agent_key is not None:
            relevant_keys |= set(ontology.get_parents(*agent_key))
    relevant_keys &= all_keys_for_role
    return relevant_keys

