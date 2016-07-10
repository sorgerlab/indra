import sys
import logging
import pygraphviz as pgv
import itertools
from copy import copy, deepcopy
from indra.statements import *
from indra.databases import uniprot_client

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
    def __init__(self, hierarchies, stmts=None):
        self.hierarchies = hierarchies
        if stmts:
            self.stmts = deepcopy(stmts)
        else:
            self.stmts = []
        self.unique_stmts = []
        self.related_stmts = []

    def add_statements(self, stmts):
        """Add to the current list of statements.

        Parameters
        ----------
        stmts : list of :py:class:`indra.statements.Statement`
            Statements to add to the current list.
        """
        self.stmts += deepcopy(stmts)

    def combine_duplicates(self):
        """Combine duplicates among `stmts` and save result in `unique_stmts`.

        A wrapper around the static method :py:meth:`combine_duplicate_stmts`.
        """
        self.unique_stmts = self.combine_duplicate_stmts(self.stmts)
        return self.unique_stmts

    @staticmethod
    def combine_duplicate_stmts(stmts):
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

        >>> map2k1 = Agent('MAP2K1')
        >>> mapk1 = Agent('MAPK1')
        >>> stmt1 = Phosphorylation(map2k1, mapk1, 'T', '185',
        ... evidence=[Evidence(text='evidence 1')])
        >>> stmt2 = Phosphorylation(map2k1, mapk1, 'T', '185',
        ... evidence=[Evidence(text='evidence 2')])
        >>> uniq_stmts = Preassembler.combine_duplicate_stmts([stmt1, stmt2])
        >>> uniq_stmts
        [Phosphorylation(MAP2K1(), MAPK1(), T, 185)]
        >>> sorted([e.text for e in uniq_stmts[0].evidence])
        ['evidence 1', 'evidence 2']
        """
        unique_stmts = []
        # Remove exact duplicates using a set() call, then make copies:
        st = list(deepcopy(set(stmts)))
        # Group statements according to whether they are matches (differing
        # only in their evidence).
        # Sort the statements in place by matches_key()
        st.sort(key=lambda x: x.matches_key())
        for key, duplicates in itertools.groupby(st,
                                                 key=lambda x: x.matches_key()):
            # Get the first statement and add the evidence of all subsequent
            # Statements to it
            for stmt_ix, stmt in enumerate(duplicates):
                if stmt_ix == 0:
                    first_stmt = stmt
                else:
                    first_stmt.evidence += stmt.evidence
            # This should never be None or anything else
            assert isinstance(first_stmt, Statement)
            unique_stmts.append(first_stmt)
        return unique_stmts

    def combine_related(self):
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

        1. The statements are grouped by whether they are of the same type
           (e.g., Phosphorylation) and involve the same entities (e.g.,
           BRAF and MAP2K1).
        2. The groups of statements are then further compared to see if one
           group involves superfamily entities of another group. If so, the
           statements involving the superfamily are added into the group of
           statements involving the more concrete entities to create an
           "extended group."
        3. The statements within each extended group are then compared; if one
           statement represents a refinement of the other (as defined by the
           `refinement_of()` method implemented for the Statement), then the
           more refined statement is added to the `supports` field of the more
           general statement, and the more general statement is added to the
           `supported_by` field of the more refined statement.
        4. A new flat list of statements is created that contains only those
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

        Returns
        -------
        list of :py:class:`indra.statement.Statement`
            The returned list contains Statements representing the more
            concrete/refined versions of the Statements involving particular
            entities. The attribute :py:attr:`related_stmts` is also set to
            this list.

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
        Combining ...
        >>> combined_stmts
        [Phosphorylation(BRAF(), MAP2K1(), S)]
        >>> combined_stmts[0].supported_by
        [Phosphorylation(BRAF(), MAP2K1())]
        >>> combined_stmts[0].supported_by[0].supports
        [Phosphorylation(BRAF(), MAP2K1(), S)]
        """
        # If unique_stmts is not initialized, call combine_duplicates.
        if not self.unique_stmts:
            self.combine_duplicates()
        # Group statements according to whether they have matching entities,
        # and store the resulting lists in a dict, indexed by the key defined
        # by the statement type and its entities:
        # Sort the statements in place by entities_match_key():
        unique_stmts = deepcopy(self.unique_stmts)
        unique_stmts.sort(key=lambda x: x.entities_match_key())
        groups = {grouper[0]: list(grouper[1])
                  for grouper in itertools.groupby(unique_stmts,
                                          key=lambda x: x.entities_match_key())}
        # The ext_groups dict is where we store the extended groups, those
        # statements which involve either the same entities or entities with
        # family relationships.
        ext_groups = copy(groups)
        # Set up progress bar
        # see http://stackoverflow.com/questions/3160699/python-progress-bar
        toolbar_width = 40
        sys.stdout.write("Combining related stmts: [%s]"
                         % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width+1)) # return to start of bar
        dashes_printed = 0
        comparisons = list(itertools.permutations(groups.keys(), 2))
        num_comparisons = len(comparisons)
        # We examine pairs of Statement groups, looking for "isa" relationships:
        for counter, (g1_key, g2_key) in enumerate(comparisons):
            # Update progress bar
            pct_completed = (counter + 1) / float(num_comparisons)
            total_dashes = int(pct_completed * toolbar_width)
            dashes_to_print = total_dashes - dashes_printed
            sys.stdout.write("-" * dashes_to_print)
            sys.stdout.flush()
            dashes_printed += dashes_to_print
            # Get the groups
            g1 = groups[g1_key]
            g2 = groups[g2_key]
            # If we have two groups G1 and G2, each containing Statements with
            # some number of Agent arguments, e.g.  G1_Stmt(Ag1, Ag2, Ag3) and
            # G2_Stmt(Ag4, Ag5, Ag6), we need to know whether each of the
            # corresponding Agent arguments are related for all of the
            # Statements in the two groups.  That is, each of the arguments
            # between G1 and G2 must either be an entity match (as determined
            # by the entities_match test of the respective agents) or have an
            # "isa" relationship, as determined by the HierarchyManager in use.
            # If the elements in both groups are related, then the
            # Statements in the group with the superfamily relationship
            # are added into the group with the more specific entities. However,
            # the superfamily group is not removed from the group list, as
            # it may be a superfamily supporting Statements from another group.

            # Get the first statement from each group (we've already determined
            # that all of the statements within each group have an entity
            # match, so we only need the first:
            g1_stmt = g1[0]
            g2_stmt = g2[0]
            # Check that the statements are of the same type; if not, no merge.
            if type(g1_stmt) is not type(g2_stmt):
                continue
            # If both statements are Complexes, make sure they have the same
            # number of members:
            if type(g1_stmt) is Complex and \
               len(g1_stmt.members) != len(g2_stmt.members):
                continue
            # Check that all of the agents match or have an isa relationship.
            # Because the statements are of the same type, they should have the
            # same number of agents as arguments.  First, let's keep track of
            # our checks that g1 is the "primary" group, i.e., it is the one
            # with the more refined/grounded entities. We build a list with one
            # boolean entry for each argument, where a True value indicates
            # that the arguments at that position imply that g1 is the primary
            # group.
            agent_pairs = zip(g1_stmt.agent_list(), g2_stmt.agent_list())
            g1_is_refinement = []
            for ag1, ag2 in agent_pairs:
                if ag2 is None:
                    val = True
                elif ag2 is not None and ag1 is None:
                    val = False
                else:
                    val = ag1.entity_matches(ag2) or\
                          self.hierarchies['entity'].isa(ag1.name, ag2.name)
                g1_is_refinement.append(val)
            # If g1_is_refinement is all True values, that means everything in
            # the group1 statements isa thing in the group2 statements.
            if all(g1_is_refinement):
                g1_ext_list = ext_groups[g1_key]
                ext_groups[g1_key] = g1_ext_list + g2
        # Move cursor to next line after progress bar
        print
        # At this point we have, in ext_groups, a dict of lists of Statements
        # indexed by their entity_matches key, but now the groups contain not
        # only statements with matching entities, but also entities related by
        # isa relationships. The next step is to process each group, checking
        # each statement against each other statement to determine supports and
        # supported_by relationships. This is determined by calling the
        # refinement_of method on pairs of statements.
        # Iterate over each of the extended groups:
        for ext_group in ext_groups.values():
            # Iterate over pairs of statements in the group:
            for stmt1, stmt2 in itertools.permutations(ext_group, 2):
                if stmt1.refinement_of(stmt2, self.hierarchies):
                    stmt1.supported_by.append(stmt2)
                    stmt2.supports.append(stmt1)
        # Now that the groups have been processed, we need to find the
        # non-subsumed statements, those that have no supports relationships.
        self.related_stmts = [stmt for ext_group in ext_groups.values()
                               for stmt in ext_group
                               if not stmt.supports]
        self.related_stmts = self.combine_duplicate_stmts(self.related_stmts)

        # Make sure we haven't lost any statements!
        assert len(flatten_stmts(self.unique_stmts)) == \
               len(flatten_stmts(self.related_stmts)), \
               "Statements lost after combining related"
        return self.related_stmts


def render_stmt_graph(statements, agent_style=None):
    """Render the statement hierarchy as a pygraphviz graph.

    Parameters
    ----------
    stmts : list of :py:class:`indra.statements.Statement`
        A list of top-level statements with associated supporting statements
        resulting from building a statement hierarchy with
        :py:meth:`combine_related`.
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
    Combining ...
    [Phosphorylation(BRAF(), MAP2K1(), S)]
    >>> graph = render_stmt_graph(pa.related_stmts)
    >>> graph.write('example_graph.dot') # To make the DOT file
    >>> graph.draw('example_graph.png', prog='dot') # To make an image

    Resulting graph:

    .. image:: /images/example_graph.png
        :align: center
        :alt: Example statement graph rendered by Graphviz

    """
    # Set the default agent formatting properties
    if agent_style is None:
        agent_style = {'color': 'lightgray', 'style': 'filled',
                       'fontname': 'arial'}
    # Sets to store all of the nodes and edges as we recursively process all
    # of the statements
    nodes = set([])
    edges = set([])
    # Recursive function for processing all statements
    def process_stmt(stmt):
        nodes.add(stmt)
        for sby_ix, sby_stmt in enumerate(stmt.supported_by):
            edges.add((str(stmt.matches_key()), str(sby_stmt.matches_key())))
            process_stmt(sby_stmt)
    # Process all of the top-level statements, getting the supporting statements
    # recursively
    for stmt in statements:
        process_stmt(stmt)
    # Add the nodes and edges to the graph
    graph = pgv.AGraph(name='statements', directed=True, rankdir='LR')
    for node in nodes:
        graph.add_node(str(node.matches_key()), label=str(node), **agent_style)
    graph.add_edges_from(edges)
    return graph


def flatten_stmts(stmts):
    """Return the full set of unique stms in a pre-assembled stmt graph.

    The flattened list of of statements returned by this function can be
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
    Combining ...
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


def _flatten_evidence_for_stmt(stmt):
    total_evidence = set(stmt.evidence)
    for supp_stmt in stmt.supported_by:
        child_evidence = _flatten_evidence_for_stmt(supp_stmt)
        total_evidence = total_evidence.union(child_evidence)
    return list(total_evidence)


def flatten_evidence(stmts):
    """Add evidence from *supporting* stmts to evidence for *supported* stmts.

    Parameters
    ----------
    stmts : list of :py:class:`indra.statements.Statement`
        A list of top-level statements with associated supporting statements
        resulting from building a statement hierarchy with
        :py:meth:`combine_related`.

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
    Combining ...
    [Phosphorylation(BRAF(), MAP2K1(), S)]
    >>> [e.text for e in pa.related_stmts[0].evidence]
    ['baz', 'bak']
    >>> flattened = flatten_evidence(pa.related_stmts)
    >>> sorted([e.text for e in flattened[0].evidence])
    ['bak', 'bar', 'baz', 'foo']
    """
    # Copy all of the statements--these will be the ones where we update
    # the evidence lists
    copied_stmts = deepcopy(stmts)
    for stmt in stmts:
        total_evidence = _flatten_evidence_for_stmt(stmt)
        stmt.evidence = total_evidence
    return stmts
