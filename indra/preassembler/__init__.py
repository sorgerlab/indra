import itertools
from indra.statements import *
from indra.preassembler.hierarchy_manager import HierarchyManager

hm = HierarchyManager('../preassembler/entity_hierarchy.rdf')

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
        """Merge related statements into a supports tree structure."""
        # Group statements according to whether they have matching entities:
        group_list = [list(grouper[1])
                      for grouper in itertools.groupby(stmts,
                                          key=lambda x: x.entities_match_key())]
        # Now, make pairs of groups:
        for g1, g2 in itertools.combinations(group_list, 2):

            # Now we compare groups. If we have two groups G1 and G2, each
            # containing Statements with some number of Agent arguments, e.g.
            # G1_Stmt(Ag1, Ag2, Ag3) and G2_Stmt(Ag4, Ag5, Ag6), we need to
            # know whether each of the corresponding Agent arguments are
            # related for all of the Statements in the two groups.  That is,
            # each of the arguments between G1 and G2 must either be an entity
            # match (as determined by the entities_match test of the respective
            # agents) or have an "isa" relationship, as determined by the
            # HierarchyManager in use.

            # If the elements in both groups are related, then the
            # elements in the group with the superfamily relationship
            # is added into the group with the more specific entities. However,
            # the superfamily group is not removed from the group list, as
            # it may be a superfamily supporting Statements from another group.

            # For now, we will require that merges can only occur if the
            # isa relationships are all in the same direction for all the
            # agents. For example, for the two statement groups below
            #     RAF_family -> MEK1
            #     BRAF -> MEK_family
            # would not be merged, since BRAF isa RAF_family, but MEK_family
            # is not a MEK1. In the future this restriction could be revisited.

            # Get the first statement from each group (we've already determined
            # that all of the statements within each group have an entity match
            g1_stmt = g1[0]
            g2_stmt = g2[0]
            # Check that the statements are of the same type; if not, no merge.
            if type(g1_stmt) != type(g2_stmt):
                continue
            # Check that all of the agents match. Because the statements are
            # of the same type, they should have the same number of agents
            # as arguments.
            # First, let's keep track of our checks that g1 is the "primary"
            # group, i.e., it is the one with the more refined/grounded
            # entities. We build a list with one boolean entry for each
            # argument, where a True value indicates that the arguments
            # at that position imply that g1 is the primary group.
            agent_pairs = zip(g1_stmt.agent_list(), g2_stmt.agent_list())
            g1_first_checks = \
                    [ag1.entity_matches(ag2) or hm.isa(ag1.name, ag2.name)
                     for ag1, ag2 in agent_pairs]
            # Now we do the same in reverse, checking for g2 isa g1:
            g2_first_checks = \
                    [ag1.entity_matches(ag2) or hm.isa(ag2.name, ag1.name)
                     for ag1, ag2 in agent_pairs]
            # Now we see what we have. If g1_first_checks is all True values,
            # that means everything in the group1 statements isa thing in
            # the group2 statements.
            print g1_stmt
            print g2_stmt
            if all(g1_first_checks):
                print "G1 primary, merge g2 into g1"
            elif all(g2_first_checks):
                print "G2 primary, merge g1 into g2"
            else:
                print "Elements do not match, do nothing."
                continue
            # After this, we end up with a list of groups, but now the groups
            # contain not only statements with matching entities, but also
            # entities related by isa relationships.

            # The next step is to process each group (perhaps in a subroutine),
            # checking each statement against each other statement to determine
            # supports supported_by relationships.

            # This will require, for each field in the statement, checking
            # whether that field is a refinement of the field in the other
            # statement.

            # In addition, during this process, the most refined element
            # at each position could be recorded in a list, so that a new
            # statement could be checked to see if it added new information
            # that was more specific than any previously observed.

            # If a statement appears that has info at a position that equal
            # in specificity to one previously observed (e.g., BRAF phos
            # MEK1 at S218, and then another stmt at S222, then perhaps the
            # multiple alternatives could be recorded in a list ['218', '222']
            # and the corresponding new statements could be generated.

            # TODO Refactor entities_match to use agent_list, and then
            # writing only a single entities_match_key method for
            # statement.
            # TODO 
            # TODO As these checks proceed, there should also be some way of
            # collecting the information at each 
            # 
            # Perhaps write a method for each Statement type that returns
            # a vector of True/False for each of the Statement args, possibly
            # also with a new Statement consisting of the most specific
            # synthesis of the two Statements.

            # When its done, there should be a few statements that support
            # no other statements. These should become the toplevel statements.
            # but are supported by other statements.
            # These should becomek

if __name__ == '__main__':
    #import pickle
    #with open('example_kami.pck') as f:
    #    stmts = pickle.load(f)

    """
    for s in stmts:
        print '-----------'
        s.print_supports()
    """
    # Two duplicate statements, both have strings, evidence should be merged
    # into lists
    # Two dupes, source has list, other has string, string is added to list
    # Two dupes source has string, other has list, list is merged
    # Two dupes, both lists, lists are merged
    # All of the above, check citations, evidence, annotations, etc.
    # What if cell type is different???

    # Add supported_by field (can be more than one obj)
    # Add supports field (can be more than one obj)

    # >2 duplicate statements, should end up with one top level, with multiple
    # supported_by (not a stack)?
