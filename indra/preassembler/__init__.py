from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import sys
import time
import logging
import itertools
import functools
import collections
import networkx as nx
import multiprocessing as mp
try:
    import pygraphviz as pgv
except ImportError:
    pass
from indra.util import fast_deepcopy
from indra.statements import *
from indra.statements import stmt_type as indra_stmt_type

logger = logging.getLogger(__name__)


class Preassembler(object):
    """De-duplicates statements and arranges them in a specificity hierarchy.

    Parameters
    ----------
    hierarchies : dict[:py:class:`indra.preassembler.hierarchy_manager`]
        A dictionary of hierarchies with keys such as 'entity' (hierarchy of
        entities, primarily specifying relationships between genes and their
        families) and 'modification' pointing to HierarchyManagers
    stmts : list of :py:class:`indra.statements.Statement` or None
        A set of statements to perform pre-assembly on. If None, statements
        should be added using the :py:meth:`add_statements` method.
    matches_fun : Optional[function]
        A functon which takes a Statement object as argument and
        returns a string key that is used for duplicate recognition. If
        supplied, it overrides the use of the built-in matches_key method of
        each Statement being assembled.
    refinement_fun : Optional[function]
        A function which takes two Statement objects and a hierarchies dict
        as an argument and returns True or False. If supplied, it overrides
        the built-in refinement_of method of each Statement being assembled.

    Attributes
    ----------
    stmts : list of :py:class:`indra.statements.Statement`
        Starting set of statements for preassembly.
    unique_stmts : list of :py:class:`indra.statements.Statement`
        Statements resulting from combining duplicates.
    related_stmts : list of :py:class:`indra.statements.Statement`
        Top-level statements after building the refinement hierarchy.
    hierarchies : dict[:py:class:`indra.preassembler.hierarchy_manager`]
        A dictionary of hierarchies with keys such as 'entity' and
        'modification' pointing to HierarchyManagers
    """
    def __init__(self, hierarchies, stmts=None, matches_fun=None,
                 refinement_fun=None):
        self.hierarchies = hierarchies
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

        >>> from indra.preassembler.hierarchy_manager import hierarchies
        >>> map2k1 = Agent('MAP2K1')
        >>> mapk1 = Agent('MAPK1')
        >>> stmt1 = Phosphorylation(map2k1, mapk1, 'T', '185',
        ... evidence=[Evidence(text='evidence 1')])
        >>> stmt2 = Phosphorylation(map2k1, mapk1, 'T', '185',
        ... evidence=[Evidence(text='evidence 2')])
        >>> pa = Preassembler(hierarchies)
        >>> uniq_stmts = pa.combine_duplicate_stmts([stmt1, stmt2])
        >>> uniq_stmts
        [Phosphorylation(MAP2K1(), MAPK1(), T, 185)]
        >>> sorted([e.text for e in uniq_stmts[0].evidence]) # doctest:+IGNORE_UNICODE
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

    def _get_entities(self, stmt, stmt_type, eh):
        entities = []
        for a in stmt.agent_list():
            # Entity is None: add the None to the entities list
            if a is None and stmt_type != Complex:
                entities.append(a)
                continue
            # Entity is not None, but could be ungrounded or not
            # in a family
            else:
                a_ns, a_id = a.get_grounding()
                # No grounding available--in this case, use the
                # entity_matches_key
                if a_ns is None or a_id is None:
                    entities.append(a.entity_matches_key())
                    continue
                # We have grounding, now check for a component ID
                uri = eh.get_uri(a_ns, a_id)
                # This is the component ID corresponding to the agent
                # in the entity hierarchy
                component = eh.components.get(uri)
                # If no component ID, use the entity_matches_key()
                if component is None:
                    entities.append(a.entity_matches_key())
                # Component ID, so this is in a family
                else:
                    # We turn the component ID into a string so that
                    # we can sort it along with entity_matches_keys
                    # for Complexes
                    entities.append(str(component))
        return entities

    def _get_stmt_by_group(self, stmt_type, stmts_this_type, eh):
        """Group Statements of `stmt_type` by their hierarchical relations."""
        # Dict of stmt group key tuples, indexed by their first Agent
        stmt_by_first = collections.defaultdict(lambda: [])
        # Dict of stmt group key tuples, indexed by their second Agent
        stmt_by_second = collections.defaultdict(lambda: [])
        # Dict of statements with None first, with second Agent as keys
        none_first = collections.defaultdict(lambda: [])
        # Dict of statements with None second, with first Agent as keys
        none_second = collections.defaultdict(lambda: [])
        # The dict of all statement groups, with tuples of components
        # or entity_matches_keys as keys
        stmt_by_group = collections.defaultdict(lambda: [])
        # Here we group Statements according to the hierarchy graph
        # components that their agents are part of
        for stmt_tuple in stmts_this_type:
            _, stmt = stmt_tuple
            entities = self._get_entities(stmt, stmt_type, eh)
            # At this point we have an entity list
            # If we're dealing with Complexes, sort the entities and use
            # as dict key
            if stmt_type == Complex:
                # There shouldn't be any statements of the type
                # e.g., Complex([Foo, None, Bar])
                assert None not in entities
                assert len(entities) > 0
                entities.sort()
                key = tuple(entities)
                if stmt_tuple not in stmt_by_group[key]:
                    stmt_by_group[key].append(stmt_tuple)
            elif stmt_type == Conversion:
                assert len(entities) > 0
                key = (entities[0],
                       tuple(sorted(entities[1:len(stmt.obj_from)+1])),
                       tuple(sorted(entities[-len(stmt.obj_to):])))
                if stmt_tuple not in stmt_by_group[key]:
                    stmt_by_group[key].append(stmt_tuple)
            # Now look at all other statement types
            # All other statements will have one or two entities
            elif len(entities) == 1:
                # If only one entity, we only need the one key
                # It should not be None!
                assert None not in entities
                key = tuple(entities)
                if stmt_tuple not in stmt_by_group[key]:
                    stmt_by_group[key].append(stmt_tuple)
            else:
                # Make sure we only have two entities, and they are not both
                # None
                key = tuple(entities)
                assert len(key) == 2
                assert key != (None, None)
                # First agent is None; add in the statements, indexed by
                # 2nd
                if key[0] is None and stmt_tuple not in none_first[key[1]]:
                    none_first[key[1]].append(stmt_tuple)
                # Second agent is None; add in the statements, indexed by
                # 1st
                elif key[1] is None and stmt_tuple not in none_second[key[0]]:
                    none_second[key[0]].append(stmt_tuple)
                # Neither entity is None!
                elif None not in key:
                    if stmt_tuple not in stmt_by_group[key]:
                        stmt_by_group[key].append(stmt_tuple)
                    if key not in stmt_by_first[key[0]]:
                        stmt_by_first[key[0]].append(key)
                    if key not in stmt_by_second[key[1]]:
                        stmt_by_second[key[1]].append(key)

        # When we've gotten here, we should have stmt_by_group entries, and
        # we may or may not have stmt_by_first/second dicts filled out
        # (depending on the statement type).
        if none_first:
            # Get the keys associated with stmts having a None first
            # argument
            for second_arg, stmts in none_first.items():
                # Look for any statements with this second arg
                second_arg_keys = stmt_by_second[second_arg]
                # If there are no more specific statements matching this
                # set of statements with a None first arg, then the
                # statements with the None first arg deserve to be in
                # their own group.
                if not second_arg_keys:
                    stmt_by_group[(None, second_arg)] = stmts
                # On the other hand, if there are statements with a matching
                # second arg component, we need to add the None first
                # statements to all groups with the matching second arg
                for second_arg_key in second_arg_keys:
                    stmt_by_group[second_arg_key] += stmts
        # Now do the corresponding steps for the statements with None as the
        # second argument:
        if none_second:
            for first_arg, stmts in none_second.items():
                # Look for any statements with this first arg
                first_arg_keys = stmt_by_first[first_arg]
                # If there are no more specific statements matching this
                # set of statements with a None second arg, then the
                # statements with the None second arg deserve to be in
                # their own group.
                if not first_arg_keys:
                    stmt_by_group[(first_arg, None)] = stmts
                # On the other hand, if there are statements with a matching
                # first arg component, we need to add the None second
                # statements to all groups with the matching first arg
                for first_arg_key in first_arg_keys:
                    stmt_by_group[first_arg_key] += stmts
        return stmt_by_group

    def _generate_id_maps(self, unique_stmts, poolsize=None,
                          size_cutoff=100, split_idx=None):
        """Connect statements using their refinement relationships."""
        # Check arguments relating to multiprocessing
        if poolsize is None:
            logger.debug('combine_related: poolsize not set, '
                         'not using multiprocessing.')
            use_mp = False
        elif sys.version_info[0] >= 3 and sys.version_info[1] >= 4:
            use_mp = True
            logger.info('combine_related: Python >= 3.4 detected, '
                        'using multiprocessing with poolsize %d, '
                        'size_cutoff %d' % (poolsize, size_cutoff))
        else:
            use_mp = False
            logger.info('combine_related: Python < 3.4 detected, '
                        'not using multiprocessing.')
        eh = self.hierarchies['entity']
        # Make a list of Statement types
        stmts_by_type = collections.defaultdict(lambda: [])
        for idx, stmt in enumerate(unique_stmts):
            stmts_by_type[indra_stmt_type(stmt)].append((idx, stmt))

        child_proc_groups = []
        parent_proc_groups = []
        skipped_groups = 0
        # Each Statement type can be preassembled independently
        for stmt_type, stmts_this_type in stmts_by_type.items():
            logger.info('Grouping %s (%s)' %
                        (stmt_type.__name__, len(stmts_this_type)))
            stmt_by_group = self._get_stmt_by_group(stmt_type, stmts_this_type,
                                                    eh)

            # Divide statements by group size
            # If we're not using multiprocessing, then all groups are local
            for g_name, g in stmt_by_group.items():
                if len(g) < 2:
                    skipped_groups += 1
                    continue
                if use_mp and len(g) >= size_cutoff:
                    child_proc_groups.append(g)
                else:
                    parent_proc_groups.append(g)

        # Now run preassembly!
        logger.debug("Groups: %d parent, %d worker, %d skipped." %
                     (len(parent_proc_groups), len(child_proc_groups),
                      skipped_groups))

        supports_func = functools.partial(_set_supports_stmt_pairs,
                                          hierarchies=self.hierarchies,
                                          split_idx=split_idx,
                                          check_entities_match=False,
                                          refinement_fun=self.refinement_fun)

        # Check if we are running any groups in child processes; note that if
        # use_mp is False, child_proc_groups will be empty
        if child_proc_groups:
            # Get a multiprocessing context
            ctx = mp.get_context('spawn')
            pool = ctx.Pool(poolsize)
            # Run the large groups remotely
            logger.debug("Running %d groups in child processes" %
                         len(child_proc_groups))
            res = pool.map_async(supports_func, child_proc_groups)
            workers_ready = False
        else:
            workers_ready = True

        # Run the small groups locally
        logger.debug("Running %d groups in parent process" %
                     len(parent_proc_groups))
        stmt_ix_map = [supports_func(stmt_tuples)
                       for stmt_tuples in parent_proc_groups]
        logger.debug("Done running parent process groups")

        while not workers_ready:
            logger.debug("Checking child processes")
            if res.ready():
                workers_ready = True
                logger.debug('Child process group comparisons successful? %s' %
                             res.successful())
                if not res.successful():
                    # The get method re-raises the underlying error that we can
                    # now catch and print.
                    try:
                        res.get()
                    except Exception as e:
                        raise Exception("Sorry, there was a problem with "
                                        "preassembly in the child processes: %s"
                                        % e)
                else:
                    stmt_ix_map += res.get()
                logger.debug("Closing pool...")
                pool.close()
                logger.debug("Joining pool...")
                pool.join()
                logger.debug("Pool closed and joined.")
            time.sleep(1)
        logger.debug("Done.")
        # Combine all redundant map edges
        stmt_ix_map_set = set([])
        for group_ix_map in stmt_ix_map:
            for ix_pair in group_ix_map:
                stmt_ix_map_set.add(ix_pair)
        return stmt_ix_map_set

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
        2. Statements of the same type are then grouped according to their
           Agents' entity hierarchy component identifiers. For instance,
           ERK, MAPK1 and MAPK3 are all in the same connected component in the
           entity hierarchy and therefore all Statements of the same type
           referencing these entities will be grouped. This grouping assures
           that relations are only possible within Statement groups and
           not among groups. For two Statements to be in the same group at
           this step, the Statements must be the same type and the Agents at
           each position in the Agent lists must either be in the same
           hierarchy component, or if they are not in the hierarchy, must have
           identical entity_matches_keys. Statements with None in one of the
           Agent list positions are collected separately at this stage.
        3. Statements with None at either the first or second position are
           iterated over. For a statement with a None as the first Agent,
           the second Agent is examined; then the Statement with None is
           added to all Statement groups with a corresponding component or
           entity_matches_key in the second position. The same procedure is
           performed for Statements with None at the second Agent position.
        4. The statements within each group are then compared; if one
           statement represents a refinement of the other (as defined by the
           `refinement_of()` method implemented for the Statement), then the
           more refined statement is added to the `supports` field of the more
           general statement, and the more general statement is added to the
           `supported_by` field of the more refined statement.
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

        >>> from indra.preassembler.hierarchy_manager import hierarchies
        >>> braf = Agent('BRAF')
        >>> map2k1 = Agent('MAP2K1')
        >>> st1 = Phosphorylation(braf, map2k1)
        >>> st2 = Phosphorylation(braf, map2k1, residue='S')
        >>> pa = Preassembler(hierarchies, [st1, st2])
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
        eh = self.hierarchies['entity']

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

            pos_stmt_by_group = self._get_stmt_by_group(pst, poss, eh)
            neg_stmt_by_group = self._get_stmt_by_group(nst, negs, eh)
            for key, pg in pos_stmt_by_group.items():
                ng = neg_stmt_by_group.get(key, [])
                for (_, st1), (_, st2) in itertools.product(pg, ng):
                    if st1.contradicts(st2, self.hierarchies):
                        contradicts.append((st1, st2))

        # Handle neutral Statements next
        neu_stmts = [Influence, ActiveForm]
        for stt in neu_stmts:
            stmts = stmts_by_type.get(stt, [])
            for (_, st1), (_, st2) in itertools.combinations(stmts, 2):
                if st1.contradicts(st2, self.hierarchies):
                    contradicts.append((st1, st2))

        return contradicts

    def _normalize_relations(self, ns, rank_key, rel_fun):
        def _replace_grounding(ns, entry, rank_key, rel_fun):
            rel_ents = rel_fun(ns, entry)
            if rel_ents:
                rel_ents = [e.split('#')[1] if '#' in e else e
                            for e in rel_ents]
                sorted_entries = sorted([entry] + rel_ents, key=rank_key)
                return sorted_entries[0]
            else:
                return entry

        if rank_key is None:
            rank_key = lambda x: x
        for stmt in self.stmts:
            for agent in stmt.agent_list():
                if agent is not None and ns in agent.db_refs:
                    grounding = agent.db_refs[ns]
                    if isinstance(grounding, list):
                        new_grounding = []
                        for entry, score in grounding:
                            new_grounding.append(
                                (_replace_grounding(ns, entry,
                                                    rank_key, rel_fun),
                                 score))
                        agent.db_refs[ns] = new_grounding
                    else:
                        agent.db_refs[ns] = _replace_grounding(ns, entry,
                                                               rank_key,
                                                               rel_fun)

    def normalize_equivalences(self, ns, rank_key=None):
        self._normalize_relations(ns, rank_key,
                                  self.hierarchies['entity'].get_equals)

    def normalize_opposites(self, ns, rank_key=None):
        self._normalize_relations(ns, rank_key,
                                  self.hierarchies['entity'].get_opposites)


def _set_supports_stmt_pairs(stmt_tuples, split_idx=None, hierarchies=None,
                             check_entities_match=False, refinement_fun=None):
    # This is useful when deep-debugging, but even for normal debug is too much.
    # logger.debug("Getting support pairs for %d tuples with idx %s and stmts "
    #              "%s split at %s."
    #              % (len(stmt_tuples), [idx for idx, _ in stmt_tuples],
    #                 [(s.get_hash(shallow=True), s) for _, s in stmt_tuples],
    #                 split_idx))
    #  Make the iterator by one of two methods, depending on the case
    if not refinement_fun:
        refinement_fun = default_refinement_fun
    if split_idx is None:
        stmt_pair_iter = itertools.combinations(stmt_tuples, 2)
    else:
        stmt_group_a = []
        stmt_group_b = []
        for idx, stmt in stmt_tuples:
            if idx <= split_idx:
                stmt_group_a.append((idx, stmt))
            else:
                stmt_group_b.append((idx, stmt))
        stmt_pair_iter = itertools.product(stmt_group_a, stmt_group_b)

    # Actually create the index maps.
    ix_map = []
    for stmt_tuple1, stmt_tuple2 in stmt_pair_iter:
        stmt_ix1, stmt1 = stmt_tuple1
        stmt_ix2, stmt2 = stmt_tuple2
        if check_entities_match and not stmt1.entities_match(stmt2):
            continue
        if refinement_fun(stmt1, stmt2, hierarchies):
            ix_map.append((stmt_ix1, stmt_ix2))
        elif refinement_fun(stmt2, stmt1, hierarchies):
            ix_map.append((stmt_ix2, stmt_ix1))
    return ix_map


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

    >>> from indra.preassembler.hierarchy_manager import hierarchies
    >>> braf = Agent('BRAF')
    >>> map2k1 = Agent('MAP2K1')
    >>> st1 = Phosphorylation(braf, map2k1)
    >>> st2 = Phosphorylation(braf, map2k1, residue='S')
    >>> pa = Preassembler(hierarchies, [st1, st2])
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

    >>> from indra.preassembler.hierarchy_manager import hierarchies
    >>> braf = Agent('BRAF')
    >>> map2k1 = Agent('MAP2K1')
    >>> st1 = Phosphorylation(braf, map2k1)
    >>> st2 = Phosphorylation(braf, map2k1, residue='S')
    >>> pa = Preassembler(hierarchies, [st1, st2])
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

    >>> from indra.preassembler.hierarchy_manager import hierarchies
    >>> braf = Agent('BRAF')
    >>> map2k1 = Agent('MAP2K1')
    >>> st1 = Phosphorylation(braf, map2k1,
    ... evidence=[Evidence(text='foo'), Evidence(text='bar')])
    >>> st2 = Phosphorylation(braf, map2k1, residue='S',
    ... evidence=[Evidence(text='baz'), Evidence(text='bak')])
    >>> pa = Preassembler(hierarchies, [st1, st2])
    >>> pa.combine_related() # doctest:+ELLIPSIS
    [Phosphorylation(BRAF(), MAP2K1(), S)]
    >>> [e.text for e in pa.related_stmts[0].evidence] # doctest:+IGNORE_UNICODE
    ['baz', 'bak']
    >>> flattened = flatten_evidence(pa.related_stmts)
    >>> sorted([e.text for e in flattened[0].evidence]) # doctest:+IGNORE_UNICODE
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


def default_refinement_fun(st1, st2, hierarchies):
    return st1.refinement_of(st2, hierarchies)


def default_matches_fun(st):
    return st.matches_key()
