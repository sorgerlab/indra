import pygraphviz as pgv
import itertools
from copy import copy
from indra.statements import *
from indra.databases import uniprot_client

class Preassembler(object):

    def __init__(self, entity_hierarchy, mod_hierarchy, stmts=None):
        self.entity_hierarchy = entity_hierarchy
        self.mod_hierarchy = mod_hierarchy
        if stmts:
            self.stmts = stmts
        else:
            self.stmts = []
        self.unique_stmts = []
        self.related_stmts = []

    def add_statements(self, stmts):
        """Add to the current list of statements."""
        self.stmts += stmts

    def combine_duplicates(self):
        self.unique_stmts = self.combine_duplicate_stmts(self.stmts)
        return self.unique_stmts

    @staticmethod
    def combine_duplicate_stmts(stmts):
        """Combine evidence from duplicate Statements.

        Statements are deemed to be duplicates if they have the same key
        returned by the matches_key() method of the Statement class.

        This function keeps the first instance of each set of duplicate
        statements and merges the evidence lists from all of the other
        statements.

        Returns a list of unique Statements and stores it in
        self.unique_stmts.
        """

        unique_stmts = []
        # Group statements according to whether they are matches (differing
        # only in their evidence).
        # Sort the statements in place by matches_key()
        stmts.sort(key=lambda x: x.matches_key())
        for key, duplicates in itertools.groupby(stmts,
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
        listed in the supports field of the top-level statements.

        If self.unique_stmts has not been initialized with the de-duplicated
        statements, self.combine_duplicates is called.

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
           refinement_of() method implemented for the Statement), then the more
           refined statement is added to the `supports` field of the more
           general statement, and the more general statement is added to the
           `supported_by` field of the more refined statement.
        4. A new flat list of statements is created that contains only those
           statements that have no `supports` entries (statements containing
           such entries are not eliminated, because they will be retrievable
           from the `supported_by` fields of other statements. This list
           is returned to the caller.

        .. note:: Subfamily relationships must be consistent across arguments

            For now, we will require that merges can only occur if the isa
            relationships are all in the same direction for all the agents in a
            Statement. For example, the two statement groups: RAF_family ->
            MEK1 and BRAF -> MEK_family would not be merged, since BRAF isa
            RAF_family, but MEK_family is not a MEK1. In the future this
            restriction could be revisited.

        Returns
        -------
        list of Statements.
            The returned list contains Statements representing the more
            concrete/refined versions of the Statements involving particular
            entities.

        Example
        -------
        A more general statement with no information about a Phosphorylation
        site is identified as supporting a more specific statement::

            >>> from indra.preassembler.hierarchy_manager import \
                    entity_hierarchy as eh, modification_hierarchy as mh
            >>> braf = Agent('BRAF')
            >>> map2k1 = Agent('MAP2K1')
            >>> st1 = Phosphorylation(braf, map2k1)
            >>> st2 = Phosphorylation(braf, map2k1, position='218')
            >>> pa = Preassembler(eh, mh, [st1, st2])
            >>> combined_stmts = pa.combine_related()
            >>> combined_stmts
            [Phosphorylation(BRAF(), MAP2K1(), None, 218)]
            >>> combined_stmts[0].supported_by
            [Phosphorylation(BRAF(), MAP2K1(), None, None)]
            >>> combined_stmts[0].supported_by[0].supports
            [Phosphorylation(BRAF(), MAP2K1(), None, 218)]
        """

        # If unique_stmts is not initialized, call combine_duplicates.
        if not self.unique_stmts:
            self.combine_duplicates()
        # Group statements according to whether they have matching entities,
        # and store the resulting lists in a dict, indexed by the key defined
        # by the statement type and its entities:
        # Sort the statements in place by entities_match_key():
        self.unique_stmts.sort(key=lambda x: x.entities_match_key())
        groups = {grouper[0]: list(grouper[1])
                  for grouper in itertools.groupby(self.unique_stmts,
                                          key=lambda x: x.entities_match_key())}
        # The ext_groups dict is where we store the extended groups, those
        # statements which involve either the same entities or entities with
        # family relationships.
        ext_groups = copy(groups)
        # We examine pairs of Statement groups, looking for "isa" relationships:
        for g1_key, g2_key in itertools.permutations(groups.keys(), 2):
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
                          self.entity_hierarchy.isa(ag1.name, ag2.name)
                g1_is_refinement.append(val)
            # If g1_is_refinement is all True values, that means everything in
            # the group1 statements isa thing in the group2 statements.
            if all(g1_is_refinement):
                g1_ext_list = ext_groups[g1_key]
                ext_groups[g1_key] = g1_ext_list + g2
            continue
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
                if stmt1.refinement_of(stmt2, self.entity_hierarchy,
                                       self.mod_hierarchy):
                    stmt1.supported_by.append(stmt2)
                    stmt2.supports.append(stmt1)
        # Now that the groups have been processed, we need to find the
        # non-subsumed statements, those that have no supports relationships.
        self.related_stmts = [stmt for ext_group in ext_groups.values()
                               for stmt in ext_group
                               if not stmt.supports]
        self.related_stmts = self.combine_duplicate_stmts(self.related_stmts) 
        return self.related_stmts

def render_stmt_graph(statements, agent_style=None):
    """Renders the supports/supported_by relationships of a set of statements
    and returns a pygraphviz graph.

    Example
    -------
    Pattern for getting statements and rendering as a Graphviz graph::

        bp = biopax_api.process_pc_pathsfromto(['BRAF'], ['MAP2K1'])
        bp.get_phosphorylation()

        pa = Preassembler(eh, mh, bp.statements)
        pa.combine_related()

        graph = render_stmt_graph(pa.related_stmts)
        graph.write('braf.dot')
        graph.draw('braf.pdf', prog='dot')
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

def check_statements(stmts, save_fname=None):
    """Iterates over a list of statements and runs checks on them. Then it
    returns a tuple of lists, with the first element containing statements
    that passed all checks, and the second the statements that failed the
    tests"""
    pass_stmts = []
    fail_stmts = []
    failures = []
    for stmt in stmts:
        print stmt
        failures += check_sequence(stmt)
        if failures:
            fail_stmts.append(stmt)
        else:
            pass_stmts.append(stmt)
    if save_fname:
        failure_set = set(failures)
        with open(save_fname, 'wt') as fh:
            for f in failure_set:
                fh.write('%s\t%s\t%s\n' % (f[0], f[1], f[2]))
    return (pass_stmts, fail_stmts)

def check_sequence(stmt):
    """Check whether references to
    residues and sequence positions are consistent with sequence
    information in the UniProt database"""
    failures = []
    if isinstance(stmt, Complex):
        for m in stmt.members:
            failures += check_agent_mod(m)
    elif isinstance(stmt, Modification):
        failures += check_agent_mod(stmt.sub)
        failures += check_agent_mod(stmt.enz)
        if stmt.position is not None:
            mc = ModCondition('phosphorylation', stmt.residue, stmt.position)
            failures += check_agent_mod(stmt.sub, [mc])
    elif isinstance(stmt, SelfModification):
        failures += check_agent_mod(stmt.sub)
        if stmt.position is not None:
            mc = ModCondition('phosphorylation', stmt.residue, stmt.position)
            failures += check_agent_mod(stmt.enz, [mc])
    elif isinstance(stmt, ActivityModification):
        failures += check_agent_mod(stmt.monomer)
        failures += check_agent_mod(stmt.monomer, stmt.mod)
    return failures

def check_agent_mod(agent, mods=None):
    failures = []
    # If no UniProt ID is found, we don't report a failure
    up_id = agent.db_refs.get('UP')
    if up_id is None:
        return failures

    # If the UniProt ID is a list then choose the first one.
    if not isinstance(up_id, basestring):
        up_id = up_id[0]

    if mods is not None:
        check_mods = mods
    else:
        check_mods = agent.mods

    for m in check_mods:
        if m.position is None:
            continue
        residue = m.residue
        if residue is None:
            continue
        ver = uniprot_client.verify_location(up_id, residue, m.position)
        if not ver:
            print '-> Sequence check failed; position %s on %s is not %s.' %\
                  (m.position, agent.name, residue)
            failures.append((agent.name, residue, m.position))
    return failures
