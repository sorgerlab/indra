import itertools
import os
from copy import copy
from indra.statements import *
from indra.preassembler.hierarchy_manager import HierarchyManager

entity_file_path = os.path.join(os.path.dirname(__file__),
                    'entity_hierarchy.rdf')
mod_file_path = os.path.join(os.path.dirname(__file__),
                    'modification_hierarchy.rdf')
eh = HierarchyManager(entity_file_path)
mh = HierarchyManager(mod_file_path)

class Preassembler(object):

    def __init__(self, stmts=None):
        if stmts:
            self.stmts = stmts
        else:
            self.stmts = []

    def add_statements(self, stmts):
        """Add to the current list of statements."""
        self.stmts += stmts

    def assemble(self):
        self.unique_stmts = Preassembler.combine_duplicates(self.stmts)
        self.related_stmts = Preassembler.combine_related(self.unique_stmts)

        # Next, refinement/subsumption relationships are considered.
        # Note that since a given statement could potentially support multiple
        # parent statements (e.g., RAF phosphorylates MEK), statements can't
        # be thrown away after the first relationship is found.
        # A list of groups of statements. Each entry is a group (a list of
        # statements of the same type, involving the same Agents).
        # If the two statements are of the same type and deal with the same
        # Agents (as determined by name or grounding), then we check
        # refinement relationships.

    @staticmethod
    def combine_duplicates(stmts):
        """Combine evidence from duplicate Statements.

        Statements are deemed to be duplicates if they have the same key
        returned by the matches_key() method of the Statement class.

        This function keeps the first instance of each set of duplicate
        statements and merges the evidence lists from all of the other
        statements.

        Returns a list of unique Statements.
        """
        unique_stmts = []
        # Group statements according to whether they are matches (differing
        # only in their evidence)
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

    @staticmethod
    def combine_related(stmts):
        """Connect related statements based on their refinement relationships.

        This function takes a flat list of statements and returns a modified
        flat list of statements containing only those statements which do not
        represent a refinement of other existing statements. In other words,
        the more general versions of a given statement do not appear at the
        top level, but instead are listed in the supports field of the
        top-level statements.

        The procedure for combining statements in this way involves a series
        of steps:

        1. The statements are grouped by whether they are of the same type
           (e.g., Phosphorylation) and involve the same entities (e.g.,
           BRAF and MAP2K1).
        2. The groups of statements are then further compared to see if one
           group involves superfamily entities of another group. If so, the
           statements involving the superfamily are added into group of
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

        Parameters
        ----------
        stmts : a list of Statements.
            The statements to be combined, preferably with duplicates removed.

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

            >>> braf = Agent('BRAF')
            >>> map2k1 = Agent('MAP2K1')
            >>> st1 = Phosphorylation(braf, map2k1, 'Phosphorylation', None)
            >>> st2 = Phosphorylation(braf, map2k1, 'Phosphorylation', '218')
            >>> combined_stmts = Preassembler.combine_related([st1, st2])
            >>> combined_stmts
            [Phosphorylation(BRAF, MAP2K1, Phosphorylation, 218, [])]
            >>> combined_stmts[0].supported_by
            [Phosphorylation(BRAF, MAP2K1, Phosphorylation, None, [])]
            >>> combined_stmts[0].supported_by[0].supports
            [Phosphorylation(BRAF, MAP2K1, Phosphorylation, 218, [])]
        """

        # Group statements according to whether they have matching entities,
        # and store the resulting lists in a dict, indexed by the key defined
        # by the statement type and its entities:
        groups = {grouper[0]: list(grouper[1])
                  for grouper in itertools.groupby(stmts,
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
            if type(g1_stmt) != type(g2_stmt):
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
            g1_is_refinement = \
                    [ag1.entity_matches(ag2) or eh.isa(ag1.name, ag2.name)
                     for ag1, ag2 in agent_pairs]
            # If g1_is_refinement is all True values, that means everything in
            # the group1 statements isa thing in the group2 statements.
            g1_ext_list = ext_groups[g1_key]
            g2_ext_list = ext_groups[g2_key]
            if all(g1_is_refinement):
                ext_groups[g1_key] = g1_ext_list + g2
                ext_groups[g2_key] = g2_ext_list
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
                if stmt1.refinement_of(stmt2, eh, mh):
                    stmt1.supported_by.append(stmt2)
                    stmt2.supports.append(stmt1)
        # Now that the groups have been processed, we need to find the
        # non-subsumed statements, those that have no supports relationships.
        combined_stmts = []
        for stmt_key, ext_group in ext_groups.iteritems():
            for stmt in ext_group:
                if not stmt.supports:
                    combined_stmts.append(stmt)

        return combined_stmts

